'''
Created on Nov 3, 2010

@author: mkiyer
'''
import argparse
import logging
import os

#import matplotlib.pyplot as plt
#import networkx as nx

from lib.transcript_parser import parse_gtf
from lib.gtf import GTFFeature
from lib.cNode import strand_int_to_str
from lib.isoform_graph import IsoformGraph, EXON
from lib.path_finder import get_isoforms

def write_bed(chrom, name, strand, score, exons):
    #print "EXONS TO PRINT", exons
    tx_start = min(e[0] for e in exons)
    tx_end = max(e[1] for e in exons)
    block_sizes = []
    block_starts = []
    for start, end in exons:
        block_starts.append(start - tx_start)
        block_sizes.append(end - start)        
    # write
    s = '\t'.join([chrom, 
                   str(tx_start), 
                   str(tx_end),
                   str(name),
                   str(score),
                   strand,
                   str(tx_start),
                   str(tx_start),
                   '0',
                   str(len(exons)),
                   ','.join(map(str,block_sizes)) + ',',
                   ','.join(map(str,block_starts)) + ','])
    return s

def write_gtf(chrom, gene_id, tx_id, strand, score, exons):
    # write transcript feature
    f = GTFFeature()
    f.seqid = chrom
    f.source = "AssemblyLine"
    f.feature_type = "transcript"
    f.start = min(e[0] for e in exons)
    f.end = max(e[1] for e in exons)
    f.score = score
    f.strand = strand
    f.phase = "."
    attrs = {}
    attrs['gene_id'] = gene_id
    attrs['tx_id'] = tx_id
    f.attrs = attrs
    return f

def find_consensus(gtf_file, overhang_threshold,
                   fraction_major_isoform,
                   max_paths):
    locus_id = 1
    gene_id = 1
    tss_id = 1
    tx_id = 1
    for locus_transcripts in parse_gtf(open(gtf_file), attr_defs=None):
        logging.debug("Locus with %d transcripts chrom=%s start=%d end=%d" % 
                      (len(locus_transcripts), locus_transcripts[0].chrom,
                       locus_transcripts[0].start, locus_transcripts[-1].end))
        isoform_graph = IsoformGraph.from_transcripts(locus_transcripts) 
        isoform_graph.collapse(overhang_threshold=overhang_threshold)        
        chrom = locus_transcripts[0].chrom        
        for loc_gene_id, loc_tss_id, score, path in get_isoforms(isoform_graph.G, 
                                                                 locus_transcripts,
                                                                 fraction_major_isoform,
                                                                 max_paths):
            gene_name = "L%07d|G%07d|TSS%07d|TU%07d" % (locus_id,
                                                        gene_id + loc_gene_id, 
                                                        tss_id + loc_tss_id,
                                                        tx_id)
            strand = strand_int_to_str(path[0].strand)            
            exons = [(node.start, node.end) for node in path]
            if strand == "-":
                exons.reverse()
            s = write_bed(chrom, gene_name, strand, score, exons)
            print s
            tx_id += 1
        gene_id += loc_gene_id + 1
        tss_id += loc_tss_id + 1 
        locus_id += 1

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")    
    parser = argparse.ArgumentParser()
    parser.add_argument("--overhang", type=int, dest="overhang_threshold", default=100)
    parser.add_argument("--fraction-major-isoform", dest="fraction_major_isoform", type=float, default=0.15)
    parser.add_argument("--max-paths", dest="max_paths", type=int, default=10)
    parser.add_argument("--score-attr", dest="score_attr", default="FPKM")
    parser.add_argument("filename")
    options = parser.parse_args()
    find_consensus(options.filename, 
                   options.overhang_threshold,
                   options.fraction_major_isoform,
                   options.max_paths)

if __name__ == '__main__': main()