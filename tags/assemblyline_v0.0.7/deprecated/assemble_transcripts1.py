'''
Created on Dec 2, 2011

@author: mkiyer
'''
import argparse
import logging
import sys

from assemblyline.lib.transcript_parser import parse_gtf, cufflinks_attr_defs
from assemblyline.lib.transcript import Exon, strand_int_to_str, NEG_STRAND, POS_STRAND, NO_STRAND
from assemblyline.lib.transcript_graph import create_transcript_graph
from assemblyline.lib.assembler_base import GLOBAL_LOCUS_ID
from assemblyline.lib.assembler import assemble_transcript_graph
from assemblyline.lib.gtf import GTFFeature

GLOBAL_LOCUS_ID = 1

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

def assemble_locus(transcripts, overhang_threshold, fraction_major_isoform, max_paths,
                   bed_fileh, gtf_fileh=None):
    # gather properties of locus
    locus_chrom = transcripts[0].chrom
    locus_start = transcripts[0].start
    locus_end = max(tx.end for tx in transcripts)
    logging.debug("[LOCUS] %s:%d-%d %d transcripts" % 
                  (locus_chrom, locus_start, locus_end, 
                   len(transcripts)))
    locus_id_str = "L%d" % (GLOBAL_LOCUS_ID)
    # build transcript graph
    GG, partial_path_lists = create_transcript_graph(transcripts, overhang_threshold)
    # assemble transcripts on each strand
    features = []
    for strand in (POS_STRAND, NEG_STRAND, NO_STRAND):
        G = GG[strand]
        partial_paths = partial_path_lists[strand]
        for gene_id, path_info_list in assemble_transcript_graph(G, strand, partial_paths, fraction_major_isoform, max_paths):
            logging.debug("[LOCUS][STRAND] %s assembled %d transcripts" % 
                          (strand_int_to_str(strand), len(path_info_list))) 
            gene_id_str = "G%d" % (gene_id)
            # compute total density in path and relative density of each transcript
            total_density = sum(p.density for p in path_info_list)
            if total_density == 0.0:
                total_density = 1e-8
            # create GTF features for each transcript path
            for p in path_info_list:
                tss_id_str = "TSS%d" % p.tss_id
                tx_id_str = "TU%d" % p.tx_id
                frac = p.density / total_density
                # write to GTF
                if gtf_fileh is not None:
                    features.extend(get_gtf_features(locus_chrom, strand, p.path,
                                                     locus_id=locus_id_str, 
                                                     gene_id=gene_id_str, 
                                                     tss_id=tss_id_str, 
                                                     transcript_id=tx_id_str,
                                                     density=p.density, 
                                                     frac=frac))                    
                # write to BED format
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

def run(gtf_file, overhang_threshold, fraction_major_isoform, max_paths,
        bed_output_file=None, gtf_output_file=None):
    global GLOBAL_LOCUS_ID
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
                       fraction_major_isoform, max_paths,
                       bed_fileh, gtf_fileh)
        GLOBAL_LOCUS_ID += 1
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
        args.fraction_major_isoform,
        args.max_paths,
        args.bed,
        args.gtf)

if __name__ == '__main__': main()