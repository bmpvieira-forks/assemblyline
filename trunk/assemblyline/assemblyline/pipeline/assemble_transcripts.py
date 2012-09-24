'''
Created on Dec 2, 2011

@author: mkiyer
'''
import argparse
import logging
import sys
import collections

from assemblyline.lib.bx.cluster import ClusterTree
from assemblyline.lib.gtf import GTFFeature
from assemblyline.lib.transcript_parser import parse_gtf, cufflinks_attr_defs
from assemblyline.lib.transcript import strand_int_to_str, NEG_STRAND
from assemblyline.lib.assemble.base import GLOBAL_LOCUS_ID, GLOBAL_GENE_ID, \
    GLOBAL_TSS_ID, GLOBAL_TRANSCRIPT_ID, NODE_SCORE
from assemblyline.lib.assemble.filter import filter_transcripts
from assemblyline.lib.assemble.transcript_graph import \
    create_transcript_graphs, prune_transcript_graph
from assemblyline.lib.assemble.assembler import assemble_transcript_graph

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

def write_bedgraph(fileh, chrom, G):
    for n in sorted(G.nodes()):
        if n.start < 0:
            continue
        fields = map(str, (chrom, n.start, n.end, G.node[n][NODE_SCORE]))
        print >>fileh, '\t'.join(fields)

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
    features = [f]
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
        features.append(f)
    return features

def annotate_gene_and_tss_ids(path_info_list, strand):
    global GLOBAL_GENE_ID
    global GLOBAL_TSS_ID
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
            tss_pos_id_map[tss_pos] = GLOBAL_TSS_ID
            tss_id = GLOBAL_TSS_ID
            GLOBAL_TSS_ID += 1
        else:
            tss_id = tss_pos_id_map[tss_pos]
        path_info.tss_id = tss_id
    # retrieve transcript clusters and assign gene ids
    for start, end, indexes in cluster_tree.getregions():
        gene_id = GLOBAL_GENE_ID
        GLOBAL_GENE_ID += 1
        for i in indexes:
            path_info_list[i].gene_id = gene_id

def assemble_gene(locus_chrom, locus_id_str, G, strand, partial_paths, 
                  sensitivity_threshold, kmax, fraction_major_isoform, 
                  max_paths, bed_fileh, gtf_fileh):
    global GLOBAL_TRANSCRIPT_ID
    # run assembly algorithm
    path_info_list = assemble_transcript_graph(G, strand, partial_paths,
                                               sensitivity_threshold, kmax,
                                               fraction_major_isoform, 
                                               max_paths)
    logging.debug("\tAssembled %d transcript(s)" % (len(path_info_list)))
    # determine gene ids and tss ids
    annotate_gene_and_tss_ids(path_info_list, strand)
    # bin transcripts by gene id
    gene_path_info_dict = collections.defaultdict(lambda: [])
    for p in path_info_list:
        gene_path_info_dict[p.gene_id].append(p)
    features = []
    for gene_path_info_list in gene_path_info_dict.itervalues():
        # highest scoring path is always first in list
        highest_score = max(1e-8, gene_path_info_list[0].score)
        # create GTF features for each transcript path
        for p in gene_path_info_list:
            # assign transcript id
            tx_id = GLOBAL_TRANSCRIPT_ID
            GLOBAL_TRANSCRIPT_ID += 1
            # get strings for each id
            tx_id_str = "TU%d" % tx_id
            tss_id_str = "TSS%d" % (p.tss_id)
            gene_id_str = "G%d" % (p.gene_id)
            # compute isoform fractions
            frac = p.score / highest_score
            # write to GTF
            if gtf_fileh is not None:
                features.extend(get_gtf_features(locus_chrom, strand, p.path,
                                                 locus_id=locus_id_str, 
                                                 gene_id=gene_id_str, 
                                                 tss_id=tss_id_str, 
                                                 transcript_id=tx_id_str,
                                                 score=p.score, 
                                                 frac=frac))                    
            # write to BED
            name = "%s|%s(%.1f)" % (gene_id_str, tx_id_str, p.score)
            fields = write_bed(locus_chrom, name, strand, int(round(1000.0*frac)), p.path)
            print >>bed_fileh, '\t'.join(fields)    
    return features

def assemble_locus(transcripts, min_transcript_length, min_trim_length, 
                   trim_utr_fraction, trim_intron_fraction, 
                   sensitivity_threshold, kmax, fraction_major_isoform, 
                   max_paths, bed_fileh, gtf_fileh=None, 
                   bedgraph_filehs=None):
    global GLOBAL_LOCUS_ID
    # gather properties of locus
    locus_chrom = transcripts[0].chrom
    locus_start = transcripts[0].start
    locus_end = max(tx.end for tx in transcripts)
    logging.debug("[LOCUS] %s:%d-%d %d transcripts" % 
                  (locus_chrom, locus_start, locus_end, 
                   len(transcripts)))
    locus_id_str = "L%d" % (GLOBAL_LOCUS_ID)
    GLOBAL_LOCUS_ID += 1
    # filter transcripts
    logging.debug("\tFiltering transcripts")
    transcripts = filter_transcripts(transcripts, min_transcript_length)    
    # build transcript graphs
    features = []
    for G, strand, strand_transcript_map in \
        create_transcript_graphs(transcripts):
        # output bedgraph data
        if bedgraph_filehs is not None:
            write_bedgraph(bedgraph_filehs[strand], locus_chrom, G)
        # process transcript graphs
        for Gsub, strand, partial_paths in \
            prune_transcript_graph(G, strand, strand_transcript_map,
                                   min_trim_length, 
                                   trim_utr_fraction,
                                   trim_intron_fraction):
            logging.debug("[GENE][STRAND] (%s) graph has %d nodes and "
                          "%d unique path fragments" % 
                          (strand_int_to_str(strand), len(Gsub), 
                           len(partial_paths)))
            # assemble subgraph
            features.extend(assemble_gene(locus_chrom, locus_id_str, 
                                          Gsub, strand, partial_paths, 
                                          sensitivity_threshold, kmax,
                                          fraction_major_isoform, 
                                          max_paths, bed_fileh, gtf_fileh))
    # output GTF
    if gtf_fileh is not None:
        # in-place sort so that 'transcript' features appear before 'exon'
        features.sort(key=lambda f: f.feature_type, reverse=True)
        features.sort(key=lambda f: f.start)
        # output transcripts to gtf
        for f in features:
            print >>gtf_fileh, str(f) 

def run(gtf_file, 
        scoring_mode,
        gtf_score_attr,
        min_transcript_length, 
        min_trim_length, 
        trim_utr_fraction, 
        trim_intron_fraction, 
        sensitivity_threshold, 
        kmax,
        fraction_major_isoform, 
        max_paths, 
        name, 
        create_bed=False, 
        create_gtf=False, 
        create_bedgraph=False):
    # setup BED output file
    if create_bed:
        bed_fileh = open(name + ".bed", "w")
    else:
        bed_fileh = sys.stdout
    # write BED file track description line
    track_priority = 1000
    bed_track_line = ['track name="%s"' % (name),
                      'description="%s"' % (name),
                      'priority=%d' % (track_priority),
                      'visibility=pack',
                      'useScore=1']
    track_priority += 1
    bed_track_line = ' '.join(bed_track_line)
    print >>bed_fileh, bed_track_line
    # setup GTF output file
    if create_gtf:
        gtf_fileh = open(name + ".gtf", "w")
    else:
        gtf_fileh = None
    # setup bedgraph output files
    if create_bedgraph:
        bedgraph_filehs = []
        strand_colors = ['255,0,0', '0,0,255', '0,0,0']
        strand_names = ['pos', 'neg', 'none']
        # write bedgraph file track description lines
        for strand in xrange(0,3):
            track_name = '%s_%s' % (name, strand_names[strand])
            fileh = open(track_name + ".bedgraph", 'w')
            bedgraph_filehs.append(fileh)
            track_line = ['track type=bedGraph',
                          'name="%s"' % (track_name),
                          'description="%s"' % (track_name),
                          'visibility=full',
                          'color=%s' % (strand_colors[strand]),
                          'priority=%d' % (track_priority + strand),
                          'autoScale=on',
                          'alwaysZero=on',
                          'maxHeightPixels=64:64:11']
            track_line = ' '.join(track_line)
            print >>fileh, track_line
    else:
        bedgraph_filehs = None
    # run assembler on each locus
    for locus_transcripts in parse_gtf(open(gtf_file), cufflinks_attr_defs):
        # assign scores to each transcript
        if scoring_mode == "recurrence":
            for t in locus_transcripts:
                t.score = 1.0
        elif scoring_mode == "gtf_attr":
            for t in locus_transcripts:
                t.score = t.attrs[gtf_score_attr]
        # assemble
        assemble_locus(locus_transcripts, 
                       min_transcript_length, 
                       min_trim_length,
                       trim_utr_fraction, 
                       trim_intron_fraction,  
                       sensitivity_threshold, 
                       kmax, 
                       fraction_major_isoform, 
                       max_paths, 
                       bed_fileh, 
                       gtf_fileh, 
                       bedgraph_filehs)
    # cleanup output files
    if create_bed:
        bed_fileh.close()
    if create_gtf:
        gtf_fileh.close()
    if create_bedgraph:
        for fileh in bedgraph_filehs:
            fileh.close()

def main():
    # parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("-v", "--verbose", action="store_true", 
                        dest="verbose", default=False)
    parser.add_argument("--scoring-mode", dest="scoring_mode", 
                        choices=("recurrence", "gtf_attr"), 
                        default="gtf_attr", metavar="MODE",
                        help="define how assembler weights nodes "
                        " [default=%(default)s]")
    parser.add_argument("--gtf-score-attr", dest="gtf_score_attr", 
                        default="FPKM", metavar="ATTR",
                        help="GTF attribute field containing node weight "
                        " [default=%(default)s]")
    parser.add_argument("--min-transcript-length", 
                        dest="min_transcript_length", type=int,  
                        default=100, metavar="N",
                        help="Ignore transcripts that are less than N bp "
                        "in length [default=%(default)s]")
    parser.add_argument("--min-trim-length", dest="min_trim_length", 
                        type=int, default=25, metavar="N",
                        help="Constrain UTR and intron trimming such that "
                        "ends of transcripts are >N bases "
                        "[default=%(default)s]")
    parser.add_argument("--trim-utr-fraction", dest="trim_utr_fraction", 
                        type=float, default=0.1, metavar="FRAC",
                        help="Trim UTRs when coverage is less than or equal"
                        "to FRAC fraction of the remaining portion of the "
                        "UTR [default=%(default)s]")
    parser.add_argument("--trim-intron-fraction", dest="trim_intron_fraction", 
                        type=float, default=0.25, metavar="FRAC",
                        help="Trim intronic coverage when less than or "
                        "equal to FRAC fraction of the downstream exon "
                        "[default=%(default)s]")
    parser.add_argument("--kmax", dest="kmax", 
                        type=int, default=None, metavar="k",
                        help="Constrain de Bruijn graph parameter 'k' "
                        "to improve runtime performance "
                        "[default=%(default)s]")
    parser.add_argument("--sensitivity", dest="sensitivity_threshold", 
                        type=float, default=0.9, metavar="X",
                        help="Constrain de Bruijn graph parameter 'k' "
                        "such that at least X fraction of total expression "
                        "is retained in k-mer graphs [default=%(default)s]")
    parser.add_argument("--fraction-major-isoform", 
                        dest="fraction_major_isoform", type=float, 
                        default=0.001, metavar="FRAC",
                        help="Report transcript isoforms with expression "
                        "fraction >=FRAC (0.0-1.0) relative to the major "
                        "isoform [default=%(default)s]")
    parser.add_argument("--max-paths", dest="max_paths", type=int, 
                        default=1000, metavar="N",
                        help="Maximum path finding iterations to perform "
                        "for each gene [default=%(default)s]")
    parser.add_argument("-n", "--name", dest="name", default="assemblyline",
                        help="name of assembly used for output files "
                        "[default=%(default)s]")
    parser.add_argument("--bed", action="store_true", dest="bed", 
                        default=False,
                        help="Produce BED output file [default behavior "
                        "is to write to standard output]")
    parser.add_argument("--gtf", action="store_true", dest="gtf", 
                        default=False,
                        help="Produce GTF output file [no file created "
                        "by default]")
    parser.add_argument("--bedgraph", action="store_true", dest="bedgraph", 
                        default=False,
                        help="Produce BedGraph track file(s) showing the "
                        "coverage profile across all transcript nodes [no "
                        "file created by default]")
    parser.add_argument("gtf_input_file")
    args = parser.parse_args()
    # set logging level
    if args.verbose:
        level = logging.DEBUG
    else:
        level = logging.WARNING
    logging.basicConfig(level=level,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    logging.info("AssemblyLine")
    logging.info("----------------------------------")
    logging.info("Parameters:")
    logging.info("input file:              %s" % (args.gtf_input_file))
    logging.info("scoring mode:            %s" % (args.scoring_mode))
    logging.info("gtf score attribute:     %s" % (args.gtf_score_attr))
    logging.info("assembly name:           %s" % (args.name))
    logging.info("min transcript length:   %d" % (args.min_transcript_length))
    logging.info("min trim length:         %d" % (args.min_trim_length))
    logging.info("trim utr fraction:       %f" % (args.trim_utr_fraction))
    logging.info("trim intron fraction:    %f" % (args.trim_intron_fraction))
    logging.info("assembler sensitivity:   %f" % (args.sensitivity_threshold))
    logging.info("fraction major isoform:  %f" % (args.fraction_major_isoform))
    logging.info("max paths:               %d" % (args.max_paths))
    # constrain parameters
    if args.min_transcript_length < 0:
        parser.error("min_transcript_length < 0")
    if args.min_trim_length < 0:
        parser.error("min_trim_length < 0")
    if (args.trim_utr_fraction < 0) or (args.trim_utr_fraction > 1):
        parser.error("trim_utr_fraction out of range (0.0-1.0)")
    if (args.trim_intron_fraction < 0) or (args.trim_intron_fraction > 1):
        parser.error("trim_intron_fraction out of range (0.0-1.0)")
    if (args.sensitivity_threshold < 0) or (args.sensitivity_threshold > 1):
        parser.error("sensitivity_threshold out of range (0.0-1.0)")
    if (args.fraction_major_isoform < 0) or (args.fraction_major_isoform > 1):
        parser.error("fraction_major_isoform out of range (0.0-1.0)")
    if (args.max_paths < 1):
        parser.error("max_paths <= 0")
    # start algorithm
    run(args.gtf_input_file,
        args.scoring_mode,
        args.gtf_score_attr,
        args.min_transcript_length,
        args.min_trim_length,
        args.trim_utr_fraction,
        args.trim_intron_fraction,
        args.sensitivity_threshold,
        args.kmax,
        args.fraction_major_isoform,
        args.max_paths,
        args.name,
        args.bed,
        args.gtf,
        args.bedgraph)

if __name__ == '__main__': main()