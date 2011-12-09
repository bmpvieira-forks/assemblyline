'''
Created on Dec 2, 2011

@author: mkiyer
'''
import argparse
import logging
import os

from assemblyline.lib.transcript_parser import parse_gtf, cufflinks_attr_defs
from assemblyline.lib.transcript import Exon, strand_int_to_str, NEG_STRAND, POS_STRAND
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

#def write_gtf(chrom, gene_id, tx_id, strand, score, exons):
#    # write transcript feature
#    f = GTFFeature()
#    f.seqid = chrom
#    f.source = "AssemblyLine"
#    f.feature_type = "transcript"
#    f.start = min(e[0] for e in exons)
#    f.end = max(e[1] for e in exons)
#    f.score = score
#    f.strand = strand
#    f.phase = "."
#    attrs = {}
#    attrs['gene_id'] = gene_id
#    attrs['tx_id'] = tx_id
#    f.attrs = attrs
#    return f
#
#def write_gtf():   
#    # transcript feature
#    f = GTFFeature()
#    f.seqid = self.chrom
#    f.source = 'assemblyline'
#    f.feature_type = 'transcript'
#    f.start = self.start
#    f.end = self.end
#    f.score = 1000.0
#    f.strand = strand_int_to_str(self.strand)
#    f.phase = '.'
#    f.attrs = self.attrs
#    features = [f]
#    # exon features
#    for i,e in enumerate(self.exons):
#        f = GTFFeature()
#        f.seqid = self.chrom
#        f.source = 'assemblyline'
#        f.feature_type = 'exon'
#        f.start = e.start
#        f.end = e.end
#        f.score = 1000.0
#        f.strand = strand_int_to_str(self.strand)
#        f.phase = '.'
#        f.attrs = {}
#        f.attrs["gene_id"] = self.attrs["gene_id"]
#        f.attrs["transcript_id"] = self.attrs["transcript_id"]
#        f.attrs["exon_number"] = i
#        features.append(f)
#    return features

#def collapse_contiguous_nodes(path, strand):
#    newpath = []    
#    n = Exon(path[0].start, path[0].end)
#    for i in xrange(1, len(path)):
#        if strand == NEG_STRAND and (n.start == path[i].end):
#            n.start = path[i].start
#        elif n.end == path[i].start:
#            n.end = path[i].end
#        else:
#            newpath.append(n)
#            n = Exon(path[i].start, path[i].end)
#    newpath.append(n)
#    return newpath

def assemble_locus(transcripts, overhang_threshold, fraction_major_isoform, max_paths):
    # gather properties of locus
    locus_chrom = transcripts[0].chrom
    locus_start = transcripts[0].start
    locus_end = max(tx.end for tx in transcripts)
    logging.debug("[LOCUS] %s:%d-%d %d transcripts" % 
                  (locus_chrom, locus_start, locus_end, 
                   len(transcripts)))
    # build and refine transcript graph
    GG = create_transcript_graph(transcripts, overhang_threshold)
    # assemble transcripts on each strand
    for strand, G in enumerate(GG):
        for gene_id, path_info_list in assemble_transcript_graph(G, strand, fraction_major_isoform, max_paths):
            logging.debug("[LOCUS][STRAND] %s assembled %d transcripts" % 
                          (strand_int_to_str(strand), len(path_info_list))) 
            for p in path_info_list:                
                # use locus/gene/tss/transcript id to make gene name
                gene_name = "L%07d|G%07d|TSS%07d|TU%07d" % (GLOBAL_LOCUS_ID, gene_id, p.tss_id, p.tx_id)
                fields = write_bed(locus_chrom, gene_name, strand, p.density, p.path)
                print '\t'.join(fields)
                # collapse contiguous nodes
                # path = collapse_contiguous_nodes(p.path, strand) 
                # fix path
                #if strand == NEG_STRAND:
                #    path.reverse()           
                #fields = write_bed(locus_chrom, gene_name, strand, p.density, path)
                #print '\t'.join(fields)

def run(gtf_file, overhang_threshold, fraction_major_isoform, max_paths):    
    global GLOBAL_LOCUS_ID
    for locus_transcripts in parse_gtf(open(gtf_file), cufflinks_attr_defs):
        assemble_locus(locus_transcripts, overhang_threshold, 
                       fraction_major_isoform, max_paths)
        GLOBAL_LOCUS_ID += 1

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
                        default=0.05, metavar="FRAC",
                        help="Report transcript isoforms with expression "
                        "fraction >=FRAC (0.0-1.0) relative to the major "
                        "isoform [default=%(default)s]")
    parser.add_argument("--max-paths", dest="max_paths", type=int, 
                        default=1000, metavar="N",
                        help="Maximum path finding iterations to perform "
                        "for each gene [default=%(default)s]")
    parser.add_argument("filename")
    args = parser.parse_args()
    # set logging level
    if args.verbose:
        level = logging.DEBUG
    else:
        level = logging.WARNING
    logging.basicConfig(level=level,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")    
    # start algorithm
    run(args.filename, 
        args.overhang_threshold,
        args.fraction_major_isoform,
        args.max_paths)

if __name__ == '__main__': main()