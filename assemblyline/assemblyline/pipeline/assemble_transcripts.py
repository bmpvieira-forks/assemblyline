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
from assemblyline.lib.assemble.base import GLOBAL_LOCUS_ID, GLOBAL_GENE_ID, GLOBAL_TSS_ID, GLOBAL_TRANSCRIPT_ID
from assemblyline.lib.assemble.transcript_graph import create_transcript_graph
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

def get_gtf_features(chrom, strand, exons, locus_id, gene_id, tss_id, 
                     transcript_id, density, frac):
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
    f.score = 1000.0
    f.strand = strand_str
    f.phase = '.'
    f.attrs = {'density': '%.3f' % density,
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
        f.score = 1000.0
        f.strand = strand_str
        f.phase = '.'
        f.attrs = {'exon_number': i+1}
        f.attrs.update(attr_dict)
        features.append(f)
    return features

def _annotate_gene_and_tss_ids(path_info_list, strand):
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


def assemble_locus(transcripts, overhang_threshold, kmax, 
                   fraction_major_isoform, max_paths,
                   bed_fileh, gtf_fileh=None):
    global GLOBAL_LOCUS_ID
    global GLOBAL_TRANSCRIPT_ID
    # gather properties of locus
    locus_chrom = transcripts[0].chrom
    locus_start = transcripts[0].start
    locus_end = max(tx.end for tx in transcripts)
    logging.debug("[LOCUS] %s:%d-%d %d transcripts" % 
                  (locus_chrom, locus_start, locus_end, 
                   len(transcripts)))
    locus_id_str = "L%d" % (GLOBAL_LOCUS_ID)
    GLOBAL_LOCUS_ID += 1
    # build transcript graph
    logging.debug("\tCreating transcript graph")
    transcript_graphs = create_transcript_graph(transcripts, overhang_threshold)
    # assemble transcripts on each strand
    features = []
    for G, partial_paths, strand in transcript_graphs:
        logging.debug("[LOCUS][STRAND] (%s) graph has %d nodes and %d partial paths" % 
                      (strand_int_to_str(strand), len(G), len(partial_paths)))
        path_info_list = assemble_transcript_graph(G, strand, partial_paths, 
                                                   kmax, fraction_major_isoform, 
                                                   max_paths)
        logging.debug("\tAssembled %d transcript(s)" % (len(path_info_list)))
        # determine gene ids and tss ids
        _annotate_gene_and_tss_ids(path_info_list, strand)
        # bin transcripts by gene id
        gene_path_info_dict = collections.defaultdict(lambda: [])
        for p in path_info_list:
            gene_path_info_dict[p.gene_id].append(p)
        for gene_path_info_list in gene_path_info_dict.itervalues():
            # highest density path is always first in list
            highest_density = max(1e-8, gene_path_info_list[0].density)
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
                frac = p.density / highest_density
                # write to GTF
                if gtf_fileh is not None:
                    features.extend(get_gtf_features(locus_chrom, strand, p.path,
                                                     locus_id=locus_id_str, 
                                                     gene_id=gene_id_str, 
                                                     tss_id=tss_id_str, 
                                                     transcript_id=tx_id_str,
                                                     density=p.density, 
                                                     frac=frac))                    
                # write to BED
                name = "%s|%s(%.1f)" % (gene_id_str, tx_id_str, p.density)
                fields = write_bed(locus_chrom, name, strand, int(round(1000.0*frac)), p.path)
                print >>bed_fileh, '\t'.join(fields)
    # output GTF
    if gtf_fileh is not None:
        # in-place sort so that 'transcript' features appear before 'exon'
        features.sort(key=lambda f: f.feature_type, reverse=True)
        features.sort(key=lambda f: f.start)
        # output transcripts to gtf
        for f in features:
            print >>gtf_fileh, str(f) 

def run(gtf_file, overhang_threshold, kmax, fraction_major_isoform, max_paths,
        bed_output_file=None, gtf_output_file=None):
    # setup output files
    if bed_output_file is not None:
        bed_fileh = open(bed_output_file, "w")
    else:
        bed_fileh = sys.stdout
    if gtf_output_file is not None:
        gtf_fileh = open(gtf_output_file, "w")
    else:
        gtf_fileh = None
    # write BED file track description line
    track_line = ('track name="AssemblyLine" '
                  'description="AssemblyLine Transcripts" '
                  'visibility=pack useScore=1')
    print >>bed_fileh, track_line
    # run assembler on each locus
    for locus_transcripts in parse_gtf(open(gtf_file), cufflinks_attr_defs):
        assemble_locus(locus_transcripts, overhang_threshold, 
                       kmax, fraction_major_isoform, max_paths,
                       bed_fileh, gtf_fileh)
    # cleanup output files
    if bed_output_file is not None:
        bed_fileh.close()
    if gtf_output_file is not None:
        gtf_fileh.close()

def main():
    # parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("-v", action="store_true", dest="verbose", default=False)
    parser.add_argument("--overhang", dest="overhang_threshold", type=int,  
                        default=100, metavar="N",
                        help="Trim ends of transcripts that extend into "
                        "introns by less than or equal to N bases "
                        "[default=%(default)s]")          
    parser.add_argument("--kmax", dest="kmax", type=int, 
                        default=3, metavar="N",
                        help="Max length of partial paths stored during assembly "
                        "and path finding steps [default=%(default)s]")
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
    parser.add_argument("--bed", dest="bed", default=None,
                        help="Produce BED output file [default behavior "
                        "is to write to standard output]")
    parser.add_argument("--gtf", dest="gtf", default=None,
                        help="Produce GTF output file [no file created "
                        "by default]")
    parser.add_argument("gtf_input_file")
    args = parser.parse_args()
    # check arguments
    if args.kmax < 2:
        parser.error("'kmax' must be set >= 2")
    # set logging level
    if args.verbose:
        level = logging.DEBUG
    else:
        level = logging.WARNING
    logging.basicConfig(level=level,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")    
    # start algorithm
    run(args.gtf_input_file, 
        args.overhang_threshold,
        args.kmax,
        args.fraction_major_isoform,
        args.max_paths,
        args.bed,
        args.gtf)

if __name__ == '__main__': main()