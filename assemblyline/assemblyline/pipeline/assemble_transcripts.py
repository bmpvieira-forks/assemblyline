'''
Created on Dec 2, 2011

@author: mkiyer

AssemblyLine: transcriptome meta-assembly from RNA-Seq

Copyright (C) 2012 Matthew Iyer

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
'''
import os
import sys
import logging
import argparse
import collections
import subprocess
import shutil
from multiprocessing import Process, JoinableQueue, Value, Lock

import assemblyline
from assemblyline.lib.bx.cluster import ClusterTree
from assemblyline.lib.base import float_check_nan, GTFAttr
from assemblyline.lib.gtf import GTFFeature
from assemblyline.lib.gtf import parse_loci, merge_sort_gtf_files
from assemblyline.lib.transcript import transcripts_from_gtf_lines, \
    strand_int_to_str, NEG_STRAND

from assemblyline.lib.assemble.base import NODE_SCORE
from assemblyline.lib.assemble.filter import filter_transcripts
from assemblyline.lib.assemble.transcript_graph import create_transcript_graphs
from assemblyline.lib.assemble.assembler import assemble_transcript_graph

class LockValue(object):
    def __init__(self, initval=0):
        self.val = Value('L', initval)
        self.lock = Lock()
    def next(self):
        with self.lock:
            cur_val = self.val.value
            self.val.value += 1
            return cur_val

SCORING_MODES = ("unweighted", "gtf_attr")
STRAND_NAMES = ('pos', 'neg', 'none')
STRAND_COLORS = ('255,0,0', '0,0,255', '0,0,0')

class RunConfig(object):
    def __init__(self):
        self.gtf_input_file = None
        self.verbose = False
        self.num_processors = 1
        self.scoring_mode = "gtf_attr"
        self.gtf_score_attr = GTFAttr.PCTRANK
        self.min_transcript_length = 250
        self.min_trim_length = 25
        self.trim_utr_fraction = 0.1
        self.trim_intron_fraction = 0.25
        self.guided = False
        self.kmax = 0
        self.ksensitivity = 0.90
        self.fraction_major_isoform = 0.01
        self.max_paths = 1000
        self.output_dir = "assembly"
        self.create_gtf = True
        self.create_bed = False
        self.create_bedgraph = False
    
    def parse_args(self):
        parser = argparse.ArgumentParser()
        parser.add_argument("-v", "--verbose", action="store_true", 
                            dest="verbose", default=self.verbose)
        parser.add_argument("-p", type=int, dest="num_processors", 
                            default=self.num_processors,
                            help="Number of processes to run in parallel "
                            "[default=%(default)s]")
        parser.add_argument("--scoring-mode", dest="scoring_mode", 
                            choices=SCORING_MODES,
                            default=self.scoring_mode, metavar="MODE",
                            help="define how assembler weights nodes "
                            " [default=%(default)s]")
        parser.add_argument("--gtf-score-attr", dest="gtf_score_attr", 
                            default=self.gtf_score_attr, metavar="ATTR",
                            help="GTF attribute field containing node weight "
                            " [default=%(default)s]")
        parser.add_argument("--min-transcript-length", 
                            dest="min_transcript_length", type=int,  
                            default=self.min_transcript_length, 
                            metavar="N",
                            help="Ignore transcripts that are less than "
                            "N bp in length [default=%(default)s]")
        parser.add_argument("--min-trim-length", dest="min_trim_length", 
                            type=int, 
                            default=self.min_trim_length, 
                            metavar="N",
                            help="Constrain UTR and intron trimming such that "
                            "ends of transcripts are >N bases "
                            "[default=%(default)s]")
        parser.add_argument("--trim-utr-fraction", 
                            dest="trim_utr_fraction", 
                            type=float, 
                            default=self.trim_utr_fraction, metavar="FRAC",
                            help="Trim UTRs when coverage is less than or equal"
                            "to FRAC fraction of the remaining portion of the "
                            "UTR [default=%(default)s]")
        parser.add_argument("--trim-intron-fraction", 
                            dest="trim_intron_fraction", 
                            type=float, 
                            default=self.trim_intron_fraction, 
                            metavar="FRAC",
                            help="Trim intronic coverage when less than or "
                            "equal to FRAC fraction of the downstream exon "
                            "[default=%(default)s]")
        grp = parser.add_argument_group("Assembly options")
        grp.add_argument("--guided", dest="guided", action="store_true",
                         default=self.guided,
                         help="Use reference transcripts as a template "
                         "during assembly (default: not set)")
        grp.add_argument("--kmax", dest="kmax", 
                         type=int, default=self.kmax, metavar="k",
                         help="Set maximum complexity of transcript "
                         "overlap graph to improve runtime performance. "
                         "Setting to zero places no limit on complexity "
                         "and may impact runtime performance on huge "
                         "datasets [default=%(default)s]")
        grp.add_argument("--ksensitivity", dest="ksensitivity", 
                         type=float, default=self.ksensitivity, 
                         metavar="X",
                         help="Optimize complexity of transcript graph "
                         "construction under the constraint that no more "
                         "than X fraction of total coverage is retained. "
                         "Setting to zero disables the optimization process "
                         "and sets 'k' equal to '--kmax' [default=%(default)s]")
        grp.add_argument("--fraction-major-isoform", 
                         dest="fraction_major_isoform", type=float, 
                         default=self.fraction_major_isoform, 
                         metavar="FRAC",
                         help="Report transcript isoforms with expression "
                         "fraction >=FRAC (0.0-1.0) relative to the major "
                         "isoform [default=%(default)s]")
        grp.add_argument("--max-paths", dest="max_paths", type=int, 
                         default=self.max_paths, metavar="N",
                         help="Maximum path finding iterations to perform "
                         "for each gene [default=%(default)s]")
        grp = parser.add_argument_group("Output options")
        grp.add_argument("-o", "--output-dir", dest="output_dir", 
                         default=self.output_dir,
                         help="directory where output files will be stored "
                         "(created if it does not exist) "
                         "[default=%(default)s]")
        grp.add_argument("--gtf", action="store_true", dest="create_gtf", 
                         default=self.create_gtf,
                         help="Produce GTF output file "
                         "[default=%(default)s]")
        grp.add_argument("--bed", action="store_true", dest="create_bed", 
                         default=self.create_bed,
                         help="Produce BED output file "
                         "[default=%(default)s]")
        grp.add_argument("--bedgraph", action="store_true", 
                         dest="create_bedgraph", 
                         default=self.create_bedgraph,
                         help="Produce bedgraph output files "
                         "[default=%(default)s]")
        parser.add_argument("gtf_input_file")
        # parse command line
        args = parser.parse_args()
        # constrain parameters
        if args.min_transcript_length < 0:
            parser.error("min_transcript_length < 0")
        if args.min_trim_length < 0:
            parser.error("min_trim_length < 0")
        if (args.trim_utr_fraction < 0) or (args.trim_utr_fraction > 1):
            parser.error("trim_utr_fraction out of range (0.0-1.0)")
        if (args.trim_intron_fraction < 0) or (args.trim_intron_fraction > 1):
            parser.error("trim_intron_fraction out of range (0.0-1.0)")
        if (args.ksensitivity < 0) or (args.ksensitivity > 1):
            parser.error("ksensitivity out of range (0.0-1.0)")
        if (args.ksensitivity < 1e-8) and (args.kmax == 0):
            parser.error("when ksensitivity set to zero please specify 'kmax' >= 1")
        if (args.kmax < 0):
            parser.error("kmax must be >= 0")
        if (args.fraction_major_isoform < 0) or (args.fraction_major_isoform > 1):
            parser.error("fraction_major_isoform out of range (0.0-1.0)")
        if (args.max_paths < 1):
            parser.error("max_paths <= 0")
        # update config attributes
        self.verbose = args.verbose
        self.num_processors = args.num_processors
        self.scoring_mode = args.scoring_mode
        self.gtf_score_attr = args.gtf_score_attr
        self.min_transcript_length = args.min_transcript_length
        self.min_trim_length = args.min_trim_length
        self.gtf_input_file = args.gtf_input_file
        self.trim_utr_fraction = args.trim_utr_fraction
        self.trim_intron_fraction = args.trim_intron_fraction
        self.guided = args.guided
        self.kmax = args.kmax
        self.ksensitivity = args.ksensitivity
        self.fraction_major_isoform = args.fraction_major_isoform
        self.max_paths = args.max_paths
        self.output_dir = args.output_dir
        self.create_gtf = args.create_gtf
        self.create_bed = args.create_bed
        self.create_bedgraph = args.create_bedgraph
    
    def log(self, logging_func=logging.info):
        logging.info("AssemblyLine version %s" % (assemblyline.__version__))
        logging.info("----------------------------------")
        logging.info("input file:              %s" % (self.gtf_input_file))
        logging.info("scoring mode:            %s" % (self.scoring_mode))
        logging.info("gtf score attribute:     %s" % (self.gtf_score_attr))
        logging.info("min transcript length:   %d" % (self.min_transcript_length))
        logging.info("min trim length:         %d" % (self.min_trim_length))
        logging.info("trim utr fraction:       %f" % (self.trim_utr_fraction))
        logging.info("trim intron fraction:    %f" % (self.trim_intron_fraction))
        logging.info("guided:                  %s" % (self.guided))
        logging.info("kmax:                    %d" % (self.kmax))
        logging.info("ksensitivity:            %f" % (self.ksensitivity))
        logging.info("fraction major isoform:  %f" % (self.fraction_major_isoform))
        logging.info("max paths:               %d" % (self.max_paths))
        logging.info("output directory:        %s" % (self.output_dir))
        logging.info("bed:                     %s" % str(self.create_bed))
        logging.info("bedgraph                 %s" % str(self.create_bedgraph))
        logging.info("gtf:                     %s" % str(self.create_gtf))
        logging.info("verbose:                 %s" % str(self.verbose))
        logging.info("num_processors:          %d" % (self.num_processors))        
        logging.info("----------------------------------")

def get_gtf_features(chrom, strand, exons, locus_id, gene_id, tss_id, 
                     transcript_id, score, frac):
    tx_start = exons[0].start
    tx_end = exons[-1].end
    strand_str = strand_int_to_str(strand)
    attr_dict = {'locus_id': locus_id,
                 'gene_id': gene_id,
                 'tss_id': tss_id,
                 'transcript_id': transcript_id}
    f = GTFFeature()
    f.seqid = chrom
    f.source = 'assemblyline'
    f.feature_type = 'transcript'
    f.start = tx_start
    f.end = tx_end
    f.score = 1000.0 * int(round(frac))
    f.strand = strand_str
    f.phase = '.'
    f.attrs = {'score': '%.3f' % score,
               'frac': '%.3f' % frac}
    f.attrs.update(attr_dict)
    yield f
    for i,e in enumerate(exons):
        f = GTFFeature()
        f.seqid = chrom
        f.source = 'assemblyline'
        f.feature_type = 'exon'
        f.start = e.start
        f.end = e.end
        f.score = int(round(frac))
        f.strand = strand_str
        f.phase = '.'
        f.attrs = {'exon_number': i+1}
        f.attrs.update(attr_dict)
        yield f
        
def write_bed(chrom, name, strand, score, exons):
    assert all(exons[0].start < x.start for x in exons[1:])
    assert all(exons[-1].end > x.end for x in exons[:-1])
    tx_start = exons[0].start
    tx_end = exons[-1].end    
    block_sizes = []
    block_starts = []
    for e in exons:
        block_starts.append(e.start - tx_start)
        block_sizes.append(e.end - e.start)        
    # make bed fields
    fields = [chrom, 
              str(tx_start), 
              str(tx_end),
              str(name),
              str(score),
              strand_int_to_str(strand),
              str(tx_start),
              str(tx_start),
              '0',
              str(len(exons)),
              ','.join(map(str,block_sizes)) + ',',
              ','.join(map(str,block_starts)) + ',']
    return fields

def sort_bed(filename, output_file, tmp_dir=None):
    args = ["sort"]
    if tmp_dir is not None:
        args.extend(["-T", tmp_dir])
    args.extend(["-k1,1", "-k2,2n", filename])
    myenv = os.environ.copy()
    myenv["LC_ALL"] = "C"
    return subprocess.call(args, stdout=open(output_file, "w"), env=myenv)

def merge_sort_files(filenames, output_file, sort_func=sort_bed, 
                     tmp_dir=None):
    tmp_file = os.path.splitext(output_file)[0] + ".unsorted.gtf"
    outfh = open(tmp_file, "w")
    for filename in filenames:
        shutil.copyfileobj(open(filename), outfh)
    outfh.close()
    sort_bed(tmp_file, output_file, tmp_dir)
    os.remove(tmp_file)

def annotate_gene_and_tss_ids(path_info_list, strand,
                              gene_id_value_obj,
                              tss_id_value_obj):
    # cluster paths to determine gene ids
    cluster_tree = ClusterTree(0,1)
    # map tss positions to unique ids
    tss_pos_id_map = {}
    for i,path_info in enumerate(path_info_list):
        start = path_info.path[0].start
        end = path_info.path[-1].end
        # cluster transcript coordinates
        cluster_tree.insert(start, end, i)
        # map TSS positions to IDs
        tss_pos = end if strand == NEG_STRAND else start
        if tss_pos not in tss_pos_id_map:
            tss_id = tss_id_value_obj.next()
            tss_pos_id_map[tss_pos] = tss_id
        else:
            tss_id = tss_pos_id_map[tss_pos]
        path_info.tss_id = tss_id
    # retrieve transcript clusters and assign gene ids
    for start, end, indexes in cluster_tree.getregions():
        gene_id = gene_id_value_obj.next()
        for i in indexes:
            path_info_list[i].gene_id = gene_id

def assemble_gene(locus_chrom, locus_id_str, 
                  gene_id_value_obj, tss_id_value_obj, t_id_value_obj,
                  G, strand, partial_paths, 
                  config, gtf_fileh, bed_fileh):
    # run assembly algorithm
    path_info_list = assemble_transcript_graph(G, strand, partial_paths,
                                               config.kmax,
                                               config.ksensitivity,
                                               config.fraction_major_isoform,
                                               config.max_paths)
    logging.debug("\tAssembled %d transcript(s)" % (len(path_info_list)))
    # determine gene ids and tss ids
    annotate_gene_and_tss_ids(path_info_list, strand,
                              gene_id_value_obj,
                              tss_id_value_obj)
    # bin transcripts by gene id
    gene_path_info_dict = collections.defaultdict(lambda: [])
    for p in path_info_list:
        gene_path_info_dict[p.gene_id].append(p)
    for gene_path_info_list in gene_path_info_dict.itervalues():
        # highest scoring path is always first in list
        highest_score = max(1e-8, gene_path_info_list[0].score)
        # create GTF features for each transcript path
        for p in gene_path_info_list:
            # assign transcript id
            t_id = t_id_value_obj.next()
            # get strings for each id
            t_id_str = "TU%d" % t_id
            tss_id_str = "TSS%d" % (p.tss_id)
            gene_id_str = "G%d" % (p.gene_id)
            # compute isoform fractions
            frac = p.score / highest_score
            # write to GTF
            if config.create_gtf:
                for f in get_gtf_features(locus_chrom, strand, p.path,
                                          locus_id=locus_id_str, 
                                          gene_id=gene_id_str, 
                                          tss_id=tss_id_str, 
                                          transcript_id=t_id_str,
                                          score=p.score, 
                                          frac=frac):
                    print >>gtf_fileh, str(f)
            # write to BED
            if config.create_bed:
                name = "%s|%s(%.1f)" % (gene_id_str, t_id_str, p.score)
                fields = write_bed(locus_chrom, name, strand, 
                                   int(round(1000.0*frac)), p.path)
                print >>bed_fileh, '\t'.join(fields)    

def assemble_locus(transcripts,
                   locus_id_value_obj,
                   gene_id_value_obj,
                   tss_id_value_obj,
                   t_id_value_obj,
                   config,
                   gtf_fileh,
                   bed_fileh,
                   bedgraph_filehs):
    # gather properties of locus
    locus_chrom = transcripts[0].chrom
    locus_start = transcripts[0].start
    locus_end = max(tx.end for tx in transcripts)
    logging.debug("[LOCUS] %s:%d-%d %d transcripts" % 
                  (locus_chrom, locus_start, locus_end, 
                   len(transcripts)))
    locus_id_str = "L%d" % (locus_id_value_obj.next())
    # filter transcripts
    logging.debug("\tFiltering transcripts")
    transcripts = filter_transcripts(transcripts, 
                                     config.min_transcript_length,
                                     config.guided)
    # build transcript graphs
    transcript_graphs = \
        create_transcript_graphs(locus_chrom, transcripts, 
                                 min_trim_length=config.min_trim_length, 
                                 trim_utr_fraction=config.trim_utr_fraction,
                                 trim_intron_fraction=config.trim_intron_fraction,
                                 create_bedgraph=config.create_bedgraph,
                                 bedgraph_filehs=bedgraph_filehs)    
    for tg in transcript_graphs:
        logging.debug("Subgraph %s:%d-%d(%s) %d nodes %d paths" %
                       (locus_chrom, locus_start, locus_end,
                        strand_int_to_str(tg.strand), len(tg.Gsub),
                        len(tg.partial_paths)))
        # assemble subgraph
        assemble_gene(locus_chrom, locus_id_str, 
                      gene_id_value_obj,
                      tss_id_value_obj,
                      t_id_value_obj,
                      tg.Gsub, tg.strand, tg.partial_paths, 
                      config,
                      gtf_fileh,
                      bed_fileh)

def assembly_worker(input_queue, 
                    locus_id_value_obj,
                    gene_id_value_obj,
                    tss_id_value_obj,
                    t_id_value_obj,
                    worker_prefix,
                    config):
    # setup output files
    gtf_fileh = None
    bed_fileh = None
    bedgraph_filehs = [None, None, None]
    if config.create_gtf:
        gtf_fileh = open(worker_prefix + ".gtf", "w")
    if config.create_bed:
        bed_fileh = open(worker_prefix + ".bed", "w")
    if config.create_bedgraph:
        for strand in xrange(0,3):
            filename = '%s_%s.bedgraph' % (worker_prefix, 
                                           STRAND_NAMES[strand])
            fileh = open(filename, 'w')
            bedgraph_filehs[strand] = fileh
    # process input
    while True:
        lines = input_queue.get()
        if len(lines) == 0:
            break
        transcripts = transcripts_from_gtf_lines(lines)
        # conserve memory
        del lines
        # assign scores to each transcript
        for t in transcripts:
            if config.scoring_mode == "unweighted":
                t.score = 1.0
            elif config.scoring_mode == "gtf_attr":
                score = t.attrs.get(config.gtf_score_attr, '0')
                t.score = float_check_nan(score)
        # assemble
        assemble_locus(transcripts,
                       locus_id_value_obj,
                       gene_id_value_obj,
                       tss_id_value_obj,
                       t_id_value_obj,                       
                       config,
                       gtf_fileh,
                       bed_fileh,
                       bedgraph_filehs)
        input_queue.task_done()
    # cleanup output files
    if config.create_bed:
        bed_fileh.close()
    if config.create_gtf:
        gtf_fileh.close()
    if config.create_bedgraph:
        for fileh in bedgraph_filehs:
            fileh.close()
    input_queue.task_done()

def run_parallel(config):
    """
    runs assembly in parallel and merges output from child processes 

    config: RunConfig object
    """
    # create temp directory
    tmp_dir = os.path.join(config.output_dir, "tmp")
    if not os.path.exists(tmp_dir):
        logging.debug("Creating tmp directory '%s'" % (tmp_dir))
        os.makedirs(tmp_dir)
    # create queue
    input_queue = JoinableQueue(maxsize=config.num_processors*3)
    # shared memory values
    locus_id_value_obj = LockValue(1)
    gene_id_value_obj = LockValue(1)
    tss_id_value_obj = LockValue(1)
    t_id_value_obj = LockValue(1)
    # start worker processes
    procs = []
    worker_prefixes = []
    for i in xrange(config.num_processors):
        worker_prefix = os.path.join(tmp_dir, "worker%03d" % (i))
        worker_prefixes.append(worker_prefix)
        args = (input_queue, 
                locus_id_value_obj,
                gene_id_value_obj,
                tss_id_value_obj,
                t_id_value_obj,
                worker_prefix,
                config)
        p = Process(target=assembly_worker, args=args)
        p.daemon = True
        p.start()
        procs.append(p)
    # parse gtf file                
    for lines in parse_loci(open(config.gtf_input_file)):
        input_queue.put(lines)
        # conserve memory
        del lines
    # stop workers
    for p in procs:
        input_queue.put([])
    # close queue
    input_queue.join()
    input_queue.close()
    # join worker processes
    for p in procs:
        p.join()
    # merge gtf files
    if config.create_gtf:
        logging.info("Merging %d worker GTF files" % 
                     (config.num_processors))
        worker_gtf_files = [prefix + ".gtf" for prefix in worker_prefixes]
        output_gtf_file = os.path.join(config.output_dir, "assembly.gtf")
        merge_sort_gtf_files(worker_gtf_files, output_gtf_file, 
                             tmp_dir=tmp_dir)
        # remove worker gtf files
        for filename in worker_gtf_files:
            if os.path.exists(filename):
                os.remove(filename)
    # merge bed files
    if config.create_bed:
        logging.info("Merging %d worker BED files" % 
                     (config.num_processors))
        worker_bed_files = [p + ".bed" for p in worker_prefixes]
        output_bed_file = os.path.join(config.output_dir, "assembly.bed")
        merge_sort_files(worker_bed_files, output_bed_file, 
                         sort_func=sort_bed, 
                         tmp_dir=tmp_dir)
        # write bed file track description line
        track_name = os.path.basename(config.output_dir)
        track_line = ' '.join(['track name="%s"' % (track_name),
                               'description="%s"' % (track_name),
                               'visibility=pack',
                               'useScore=1'])
        track_file = os.path.join(config.output_dir, 
                                  "assembly.bed.ucsc_track")
        fileh = open(track_file, "w")
        print >>fileh, track_line
        fileh.close()
    # merge bedgraph files
    if config.create_bedgraph:
        logging.info("Merging %d worker bedGraph files" % 
                     (config.num_processors))
        for strand in xrange(0,3):
            strand_name = STRAND_NAMES[strand]
            bgfiles = ['%s_%s.bedgraph' % (p, strand_name)
                       for p in worker_prefixes]
            output_file = os.path.join(config.output_dir, 
                                       "assembly_%s.bedgraph" % strand_name)
            merge_sort_files(bgfiles, output_file, 
                             sort_func=sort_bed, 
                             tmp_dir=tmp_dir)
            track_name = '%s_%s' % (os.path.basename(config.output_dir), 
                                    strand_name)
            track_line = ' '.join(['track type=bedGraph',
                                   'name="%s"' % (track_name),
                                   'description="%s"' % (track_name),
                                   'visibility=full',
                                   'color=%s' % (STRAND_COLORS[strand]),
                                   'autoScale=on',
                                   'alwaysZero=on',
                                   'maxHeightPixels=64:64:11'])
            track_file = os.path.join(config.output_dir, 
                                      "assembly_%s.bedgraph.ucsc_track" % strand_name)
            fileh = open(track_file, "w")
            print >>fileh, track_line
            fileh.close()
    # cleanup
    if os.path.exists(tmp_dir):
        shutil.rmtree(tmp_dir)
    logging.info("Done")
    return 0

def main():
    # create default run configuration
    config = RunConfig()
    # parse command line
    config.parse_args()
    # setup logging
    if config.verbose:
        level = logging.DEBUG
    else:
        level = logging.INFO
    logging.basicConfig(level=level,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    # log run configuration
    config.log()
    # create output directory
    config.output_dir = os.path.abspath(config.output_dir)
    if not os.path.exists(config.output_dir):
        logging.debug("Creating output directory '%s'" % (config.output_dir))
        os.makedirs(config.output_dir)
    # start algorithm
    return run_parallel(config)

if __name__ == '__main__': 
    sys.exit(main())