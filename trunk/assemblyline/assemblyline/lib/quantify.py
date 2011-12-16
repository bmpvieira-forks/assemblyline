'''
Created on Dec 15, 2011

@author: mkiyer
'''
import collections
import logging
import numpy as np
import pysam

from bx.cluster import ClusterTree
from transcript import Exon, POS_STRAND, NEG_STRAND, NO_STRAND, strand_str_to_int, strand_int_to_str
from transcript_parser import parse_gtf

CIGAR_M = 0 #match  Alignment match (can be a sequence match or mismatch)
CIGAR_I = 1 #insertion  Insertion to the reference
CIGAR_D = 2 #deletion  Deletion from the reference
CIGAR_N = 3 #skip  Skipped region from the reference
CIGAR_S = 4 #softclip  Soft clip on the read (clipped sequence present in <seq>)
CIGAR_H = 5 #hardclip  Hard clip on the read (clipped sequence NOT present in <seq>)
CIGAR_P = 6 #padding  Padding (silent deletion from the padded reference sequence)

TAG_STRAND = "XS"
TAG_NUM_HITS = "NH"

def get_genomic_intervals(read):
    intervals = []
    astart = read.pos
    aend = astart
    for op,length in read.cigar:
        if (op == CIGAR_D) or (op == CIGAR_M):
            aend += length
        elif (op == CIGAR_N):
            if aend > astart:
                intervals.append((astart, aend))
            astart = aend + length
            aend = astart
    if aend > astart:
        intervals.append((astart, aend))
    return intervals

class Gene(object):
    def __init__(self):
        self.strand = NO_STRAND
        self.exons = set()

def quantify_locus(locus_chrom, locus_start, locus_end, reads, transcripts):
    """
    uses an exon union model to quantify gene expression
    of the 'reads' against a list of 'transcripts'
    
    generator function that returns tuples:

    (gene_id, strand, start, end, num_exons, length, fragments)
    """
    # aggregate exons for each gene (exon union model)
    gene_id_map = collections.defaultdict(lambda: Gene())
    for t in transcripts:
        g = gene_id_map[t.gene_id]
        g.exons.update(t.exons)
        g.strand = t.strand
    # cluster exons together for each gene
    A = np.zeros((3,(locus_end - locus_start)), dtype=bool)
    for g in gene_id_map.itervalues():
        cluster_tree = ClusterTree(0,1)
        for i,e in enumerate(g.exons):
            cluster_tree.insert(e.start, e.end, i)
        exon_clusters = []
        for start, end, indexes in cluster_tree.getregions():
            exon_clusters.append(Exon(start,end))
            A[g.strand,(start-locus_start):(end-locus_start)] = True
        # update exons
        g.exons = exon_clusters
        del cluster_tree
    # normalize strand fracs
    strand_fracs = np.zeros((2,(locus_end-locus_start)),dtype=float)
    # split areas with genes on both pos and neg in half
    pos_and_neg = np.logical_and(A[POS_STRAND],A[NEG_STRAND])
    pos_only = np.logical_and(A[POS_STRAND], np.logical_not(A[NEG_STRAND]))
    neg_only = np.logical_and(A[NEG_STRAND], np.logical_not(A[POS_STRAND]))
    strand_fracs[POS_STRAND,pos_and_neg] = 0.5
    strand_fracs[POS_STRAND,pos_only] = 1.0
    strand_fracs[NEG_STRAND,pos_and_neg] = 0.5
    strand_fracs[NEG_STRAND,neg_only] = 1.0
    # build numpy array of reads
    mass_arr = np.zeros((3,(locus_end - locus_start)), dtype=float)
    total_frags = 0.0
    for r in reads:
        # get read tag attributes
        tag_dict = dict(r.tags)
        strand = strand_str_to_int(tag_dict.get(TAG_STRAND, "."))
        num_hits = tag_dict.get(TAG_NUM_HITS, 1)
        # compute read coverage such that FPKM will be accurate
        intervals = get_genomic_intervals(r)
        alen = sum((end-start) for (start,end) in intervals)
        cov = 1. / (2.0 * num_hits * alen)
        # add read to mass array
        for start,end in intervals:
            mass_arr[strand, (start-locus_start):(end-locus_start)] += cov
        # keep track of total frags
        total_frags += (1./(2.0*num_hits))
    # compute the fraction of strand mass at each base and
    # update pos strand frac array
    total_stranded_mass_arr = mass_arr[POS_STRAND] + mass_arr[NEG_STRAND]
    nzinds = np.flatnonzero(total_stranded_mass_arr)
    strand_fracs[POS_STRAND,nzinds] = mass_arr[POS_STRAND,nzinds] / total_stranded_mass_arr[nzinds]
    strand_fracs[NEG_STRAND,nzinds] = mass_arr[NEG_STRAND,nzinds] / total_stranded_mass_arr[nzinds]
    # reallocate unstranded mass proportionately
    mass_delta_arr = np.zeros((locus_end - locus_start), dtype=float)
    for strand in (POS_STRAND, NEG_STRAND):
        arr = strand_fracs[strand] * mass_arr[NO_STRAND]
        mass_arr[strand] += arr
        mass_delta_arr += arr
    mass_arr[NO_STRAND] -= mass_delta_arr
    # compute abundance of genes using arrays
    gene_abundances = []
    for gene_id, g in gene_id_map.iteritems():
        # sum mass on all strands over the exonic intervals of the transcript
        frags = 0.0
        length = 0
        for e in g.exons:
            frags += np.sum(mass_arr[g.strand,(e.start - locus_start):(e.end - locus_start)])
            length += (e.end - e.start)
        frags = int(round(frags,0))
        gene_abundances.append((gene_id, strand_int_to_str(g.strand), locus_chrom, g.exons[0].start, g.exons[-1].end, len(g.exons), length, frags))
    return int(round(total_frags,0)), gene_abundances

def quantify_gene_abundance(gtf_file, bam_file, output_file):
    outfh = open(output_file, "w")
    bamfh = pysam.Samfile(bam_file, "rb")
    header = ["gene_id", "strand", "chrom", "start", "end", "num_exons", "length", "frags"]
    print >>outfh, '\t'.join(header)
    total_frags = 0
    for transcripts in parse_gtf(open(gtf_file)):
        chrom = transcripts[0].chrom
        start = transcripts[0].start
        end = max(tx.end for tx in transcripts)
        #logging.debug("[LOCUS] %s:%d-%d %d transcripts" % 
        #              (chrom, start, end, len(transcripts)))
        frags, abundances = quantify_locus(chrom, start, end, 
                                           bamfh.fetch(chrom, start, end), 
                                           transcripts)
        for itm in abundances:
            print >>outfh, '\t'.join(map(str,itm))
        total_frags += frags
    print >>outfh, "# TOTAL_FRAGS=%d" % (total_frags)
    bamfh.close()

def main():
    import argparse
    # parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("gtf_file")
    parser.add_argument("bam_file")
    parser.add_argument("output_file")
    args = parser.parse_args()
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")    
    quantify_gene_abundance(args.gtf_file, args.bam_file, args.output_file)

if __name__ == '__main__': 
    main()
