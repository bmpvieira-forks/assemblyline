'''
Created on Oct 5, 2012

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
import logging
import argparse
import os
import collections
import sys
import shutil
import numpy as np
from multiprocessing import Process, JoinableQueue

import assemblyline
from assemblyline.lib.bx.intersection import Interval, IntervalTree
from assemblyline.lib.gtf import GTFAttr, parse_loci, merge_sort_gtf_files
from assemblyline.lib.librarytable import LibraryInfo
from assemblyline.lib.transcript import parse_gtf, \
    transcripts_from_gtf_lines, NO_STRAND, POS_STRAND, NEG_STRAND
from assemblyline.lib.base import LIB_COUNTS_FILE, SENSE, ANTISENSE, \
    INTRONIC, INTERGENIC, CATEGORIES, category_int_to_str, ANNOTATED, \
    CATEGORY, MEAN_SCORE, MEAN_RECURRENCE, LibCounts, CategoryInfo

from assemblyline.lib.assemble.base import STRAND_SCORE, TRANSCRIPT_IDS
from assemblyline.lib.assemble.transcript_graph import \
    create_undirected_transcript_graph, get_transcript_node_map

# graph attributes
IS_REF = 'isref'
SAMPLE_IDS = 'sids'
# default parameters
DEFAULT_ANNOTATION_FRAC_THRESHOLD = 0.9

def get_strand_score_fraction(G, nodes):
    # sum the coverage score across the transcript
    score_arr = np.zeros(3,float)
    for n in nodes:
        score_arr += G.node[n][STRAND_SCORE]
    # calculate total mass across transcript
    total_strand_score = score_arr[POS_STRAND] + score_arr[NEG_STRAND]
    # if there is "stranded" mass on any of the nodes comprising
    # this transcript, then use the proportion of fwd/rev mass
    # to redistribute
    if total_strand_score > 0:
        # proportionally assign unstranded mass based on amount of
        # plus and minus strand mass
        pos_frac = score_arr[POS_STRAND] / float(total_strand_score)
        return pos_frac
    return None

def get_annotated_length(G, nodes, strand):
    if strand != NEG_STRAND:
        sense_strand, asense_strand = POS_STRAND, NEG_STRAND
    else:
        sense_strand, asense_strand = NEG_STRAND, POS_STRAND
    sense_ann_length = 0
    asense_ann_length = 0
    total_length = 0
    for n in nodes:
        length = (n.end - n.start)
        total_length += length
        isref = G.node[n][IS_REF]
        sense_isref = isref[sense_strand]
        asense_isref = isref[asense_strand]
        if sense_isref:
            sense_ann_length += length
        if asense_isref:
            asense_ann_length += length    
    return sense_ann_length, asense_ann_length, total_length

def get_recurrence_and_score(G, nodes, strand):
    # gather recurrence and scores across each node
    total_length = 0.0
    total_score = 0.0
    total_recur = 0.0
    for n in nodes:
        nd = G.node[n]
        length = float(n.end - n.start)
        total_score += (nd[STRAND_SCORE][strand] + nd[STRAND_SCORE][NO_STRAND]) * length
        total_recur += len(nd[SAMPLE_IDS]) * length
        total_length += length
    # calculate statistics
    mean_score = total_score / total_length
    mean_recur = total_recur / total_length
    return mean_score, mean_recur

def resolve_strand(G, nodes):
    # first try to use score fraction
    pos_frac = get_strand_score_fraction(G, nodes)
    if pos_frac is not None:
        strand = POS_STRAND if pos_frac >= 0.5 else NEG_STRAND
    else:
        # get length of annotated regions
        pos_ann_bp, neg_ann_bp, total_bp = \
            get_annotated_length(G, nodes, POS_STRAND)
        # use whichever of the two strands better matches an 
        # annotated gene
        pos_ann_frac = pos_ann_bp / float(total_bp)
        neg_ann_frac = neg_ann_bp / float(total_bp)
        if pos_ann_frac >= neg_ann_frac:
            strand = POS_STRAND
        else:
            strand = NEG_STRAND
    return strand

def add_node_annotate(G, n, t, **kwargs):
    if n not in G: 
        attr_dict = {TRANSCRIPT_IDS: set(),
                     IS_REF: np.zeros(3,bool),
                     SAMPLE_IDS: set(),
                     STRAND_SCORE: np.zeros(3,float)} 
        G.add_node(n, attr_dict=attr_dict)
    t_id = t.attrs[GTFAttr.TRANSCRIPT_ID]
    is_ref = bool(int(t.attrs[GTFAttr.REF]))
    sample_id = None
    if "gtf_sample_attr" in kwargs:
        sample_id_attr = kwargs["gtf_sample_attr"]
        if sample_id_attr in t.attrs:
            sample_id = t.attrs[sample_id_attr]
    nd = G.node[n]
    nd[TRANSCRIPT_IDS].add(t_id)
    if sample_id is not None:
        nd[SAMPLE_IDS].add(sample_id)
    nd[IS_REF][t.strand] |= is_ref
    nd[STRAND_SCORE][t.strand] += t.score

def annotate_locus(transcripts, 
                   gtf_sample_attr, 
                   gtf_score_attr,
                   annotation_frac_threshold):                   
    # get reference introns and score transcripts
    ref_introns = {POS_STRAND: set(), 
                   NEG_STRAND: set()}
    ref_transcripts = {POS_STRAND: [],
                       NEG_STRAND: []}
    intron_tree = IntervalTree()
    inp_transcripts = []
    for t in transcripts:
        # set transcript score
        t.score = float(t.attrs.get(gtf_score_attr, 0.0))
        # separate ref and nonref transcripts
        is_ref = bool(int(t.attrs[GTFAttr.REF]))
        if is_ref:
            ref_transcripts[t.strand].append(t)
            for start,end in t.iterintrons():
                ref_introns[t.strand].add((start,end))
                intron_tree.insert_interval(Interval(start,end))
        else:
            inp_transcripts.append(t)
            for start,end in t.iterintrons():
                intron_tree.insert_interval(Interval(start,end))
    # create undirected transcript graph from all transcripts
    G = create_undirected_transcript_graph(transcripts, 
                                           add_node_func=add_node_annotate,
                                           gtf_sample_attr=gtf_sample_attr)
    # build a mapping from transcripts to graph nodes using the 
    # transcript id attributes of the nodes
    transcript_node_map = get_transcript_node_map(G)
    strand_transcripts_dict = {POS_STRAND: [], NEG_STRAND: []}
    for t in inp_transcripts:
        nodes = transcript_node_map[t.attrs[GTFAttr.TRANSCRIPT_ID]]
        # resolve strand of unstranded transcripts 
        if t.strand == NO_STRAND:
            strand = resolve_strand(G, nodes)
        else:
            strand = t.strand
        # check if all introns are annotated
        myintrons = set(t.iterintrons())
        introns_annotated = True
        for intron in myintrons:
            if intron not in ref_introns[strand]:
                introns_annotated = False
                break
        # get length of annotated regions
        sense_ann_bp, asense_ann_bp, total_bp = \
            get_annotated_length(G, nodes, strand)
        # calc fraction of transcript that is annotated
        sense_ann_frac = sense_ann_bp / float(total_bp)
        is_annotated = (introns_annotated and 
                        (sense_ann_frac >= annotation_frac_threshold))
        # determine transcript category
        if is_annotated:
            category = SENSE
        elif (sense_ann_bp == 0) and (asense_ann_bp == 0):
            # search for introns overlapping transcript
            found_intron = False
            for hit in intron_tree.find(t.start, t.end):
                # ignore contained introns
                if (hit.start > t.start) and (hit.end < t.end):
                    continue
                # ignore introns of this transcript
                if (hit.start,hit.end) in myintrons:
                    continue
                found_intron = True
                break
            if found_intron:
                category = INTRONIC
            else:
                category = INTERGENIC
        elif sense_ann_bp > asense_ann_bp:
            category = SENSE
        else:
            category = ANTISENSE
        # add attributes
        t.attrs[ANNOTATED] = 1 if is_annotated else 0                                    
        t.attrs[CATEGORY] = category
        # update transcript attributes and save to list
        strand_transcripts_dict[strand].append(t)
    # annotate score and recurrence for transcripts
    for strand, strand_transcripts in strand_transcripts_dict.iteritems():
        # create undirected transcript graph for transcripts on 
        # specific strand
        G = create_undirected_transcript_graph(strand_transcripts,
                                               add_node_func=add_node_annotate,
                                               gtf_sample_attr=gtf_sample_attr)
        # build a mapping from transcripts to graph nodes using the 
        # transcript id attributes of the nodes
        transcript_node_map = get_transcript_node_map(G)
        for t in strand_transcripts:
            nodes = transcript_node_map[t.attrs[GTFAttr.TRANSCRIPT_ID]]
            # calculate recurrence and score statistics
            mean_score, mean_recur = get_recurrence_and_score(G, nodes, strand)
            t.attrs[MEAN_SCORE] = mean_score
            t.attrs[MEAN_RECURRENCE] = mean_recur
    return inp_transcripts

def annotate_gtf_worker(input_queue, gtf_file, gtf_sample_attr, 
                        gtf_score_attr, 
                        annotation_frac_threshold):
    fileh = open(gtf_file, 'w')
    while True:
        lines = input_queue.get()
        if len(lines) == 0:
            break
        transcripts = transcripts_from_gtf_lines(lines)
        annotated_transcripts = annotate_locus(transcripts, 
                                               gtf_sample_attr, 
                                               gtf_score_attr,
                                               annotation_frac_threshold) 
        for t in annotated_transcripts:
            features = t.to_gtf_features()
            for f in features:
                print >>fileh, str(f) 
        input_queue.task_done()
    fileh.close()
    input_queue.task_done()

def annotate_gtf_parallel(input_gtf_file,
                          output_gtf_file, 
                          gtf_sample_attr, 
                          gtf_score_attr, 
                          annotation_frac_threshold,                          
                          num_processors, 
                          tmp_dir):
    # create queue
    input_queue = JoinableQueue(maxsize=num_processors*3)
    # start worker processes
    procs = []
    worker_gtf_files = []
    for i in xrange(num_processors):
        worker_gtf_file = os.path.join(tmp_dir, "annotate_worker%03d.gtf" % (i))
        worker_gtf_files.append(worker_gtf_file)
        args = (input_queue, worker_gtf_file, gtf_sample_attr, gtf_score_attr, 
                annotation_frac_threshold)
        p = Process(target=annotate_gtf_worker, args=args)
        p.daemon = True
        p.start()
        procs.append(p)
    fileh = open(input_gtf_file)
    for lines in parse_loci(fileh):
        input_queue.put(lines)
    fileh.close()
    # stop workers
    for p in procs:
        input_queue.put([])
    # close queue
    input_queue.join()
    input_queue.close()
    # join worker processes
    for p in procs:
        p.join()
    # merge/sort worker gtf files
    logging.info("Merging %d worker GTF files" % (num_processors))
    merge_sort_gtf_files(worker_gtf_files, output_gtf_file, tmp_dir=tmp_dir)
    # remove worker gtf files
    for filename in worker_gtf_files:
        if os.path.exists(filename):
            os.remove(filename)

def output_by_category(gtf_file, library_ids, output_dir):
    # setup output by category
    category_info_dict = {}
    for category_key, category_str in category_int_to_str.iteritems():
        cinfo = CategoryInfo.create(library_ids, category_key, 
                                    category_str, output_dir)
        cinfo.output_gtf_fh = open(cinfo.output_gtf_file, "w")
        category_info_dict[category_key] = cinfo
    # function to gather transcript attributes
    lib_counts_dict = collections.defaultdict(lambda: LibCounts())
    for locus_transcripts in parse_gtf(open(gtf_file)):
#        locus_chrom = locus_transcripts[0].chrom
#        locus_start = locus_transcripts[0].start
#        locus_end = max(t.end for t in locus_transcripts)
#        logging.debug("[LOCUS] %s:%d-%d %d transcripts" % 
#                      (locus_chrom, locus_start, locus_end, 
#                       len(locus_transcripts)))
        # write classification table files
        category_features = collections.defaultdict(lambda: [])
        for t in locus_transcripts:
            library_id = t.attrs[GTFAttr.LIBRARY_ID]
            category = int(t.attrs[CATEGORY])
            is_annotated = int(t.attrs[ANNOTATED])
            features = t.to_gtf_features()
            # keep statistics
            lib_counts = lib_counts_dict[library_id]
            lib_counts.library_id = library_id
            lib_counts.annotated_counts[is_annotated] += 1
            if is_annotated:
                # write annotated transcripts to all categories
                for category_key in CATEGORIES:
                    category_features[category_key].extend(features)
            else:
                category_features[category].extend(features)
                # update stats
                lib_counts.category_counts[category] += 1
        # sort and output category GTF features
        for category_key, features in category_features.iteritems():
            # sort so that 'transcript' features appear before 'exon'
            features.sort(key=lambda f: f.feature_type, reverse=True)
            features.sort(key=lambda f: f.start)
            # output transcripts to gtf
            output_gtf_fh = category_info_dict[category_key].output_gtf_fh
            for f in features:
                print >>output_gtf_fh, str(f) 
    # close open file handles
    for category_key, cinfo in category_info_dict.iteritems():
        cinfo.output_gtf_fh.close()
    # write library category statistics
    logging.info("Writing library transcript count statistics")
    lib_counts_file = os.path.join(output_dir, LIB_COUNTS_FILE)
    fh = open(lib_counts_file, "w")
    print >>fh, '\t'.join(LibCounts.header_fields())
    for lib_counts in lib_counts_dict.itervalues():
        fields = lib_counts.to_fields()
        print >>fh, '\t'.join(map(str, fields))
    fh.close()

def main():
    # parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("-v", "--verbose", action="store_true", 
                        dest="verbose", default=False)
    parser.add_argument("-p", type=int, dest="num_processors", default=1,
                        help="Number of processes to run in parallel "
                        "[default=%(default)s]")
    parser.add_argument("--gtf-sample-attr", dest="gtf_sample_attr", 
                        default="sample_id", metavar="ATTR",
                        help="GTF attribute field used to distinguish "
                        "independent samples in order to compute "
                        "recurrence [default=%(default)s]")
    parser.add_argument("--gtf-score-attr", dest="gtf_score_attr", 
                        default="FPKM", metavar="ATTR",
                        help="GTF attribute field containing node weight "
                        " [default=%(default)s]")
    parser.add_argument("--annotation-frac", 
                        dest="annotation_frac_threshold", 
                        type=float, 
                        default=DEFAULT_ANNOTATION_FRAC_THRESHOLD,
                        help="fraction of transcript bases that must "
                        "overlap a reference transcript in order for "
                        "to be considered 'annotated' "
                        "[default=%(default)s]")
    parser.add_argument("-o", "--output-dir", dest="output_dir", 
                        default="transcripts",
                        help="output directory [default=%(default)s]")
    parser.add_argument("gtf_file")
    parser.add_argument("library_table_file")
    args = parser.parse_args()
    # check command line parameters
    if not os.path.exists(args.gtf_file):
        parser.error("Reference GTF file %s not found" % (args.ref_gtf_file))
    if not os.path.exists(args.library_table_file):
        parser.error("Library table file %s not found" % (args.library_table_file))
    # set logging level
    if args.verbose:
        level = logging.DEBUG
    else:
        level = logging.INFO
    logging.basicConfig(level=level,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    logging.info("AssemblyLine %s" % (assemblyline.__version__))
    logging.info("----------------------------------")   
    # show parameters
    logging.info("num processors:       %d" % (args.num_processors))
    logging.info("gtf sample attribute: %s" % (args.gtf_sample_attr))
    logging.info("gtf score attribute:  %s" % (args.gtf_score_attr))
    logging.info("output directory:     %s" % (args.output_dir))
    logging.info("gtf file:             %s" % (args.gtf_file))
    logging.info("library table file:   %s" % (args.library_table_file))
    logging.info("----------------------------------")   
    if not os.path.exists(args.output_dir):
        logging.debug("Creating output directory '%s'" % (args.output_dir))
        os.makedirs(args.output_dir)
    tmp_dir = os.path.join(args.output_dir, "tmp")
    if not os.path.exists(tmp_dir):
        logging.info("Creating tmp directory '%s'" % (tmp_dir))
        os.makedirs(tmp_dir)
    # parse library table
    logging.info("Parsing library table")
    library_ids = []
    for lib in LibraryInfo.from_file(args.library_table_file):
        library_ids.append(lib.library_id)
    # function to gather transcript attributes
    logging.info("Annotating GTF file")
    annotated_gtf_file = os.path.join(args.output_dir, "annotated_transcripts.gtf")
    annotate_gtf_parallel(args.gtf_file,
                          annotated_gtf_file,
                          args.gtf_sample_attr, 
                          args.gtf_score_attr, 
                          args.annotation_frac_threshold,                          
                          args.num_processors, 
                          tmp_dir)
    # output annotated transcripts by category
    logging.info("Outputting transcripts by category")
    output_by_category(annotated_gtf_file, library_ids, args.output_dir)
    # cleanup
    if os.path.exists(tmp_dir):
        shutil.rmtree(tmp_dir)

if __name__ == "__main__":
    sys.exit(main())
