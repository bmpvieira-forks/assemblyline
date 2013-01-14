'''
Created on Mar 13, 2011

@author: mkiyer
'''
import logging
import argparse
import collections

from lib.bx.intersection import Interval, IntervalTree
from lib.bx.cluster import ClusterTree
from lib.transcript_parser import parse_gtf
from lib.base import POS_STRAND, NEG_STRAND, NO_STRAND
import lib.bed as bed

def build_interval_trees(bed_iter):
    exons = collections.defaultdict(lambda: set())
    introns = collections.defaultdict(lambda: set())
    exon_trees = collections.defaultdict(lambda: IntervalTree())
    intron_trees = collections.defaultdict(lambda: IntervalTree())
    for g in bed.BEDFeature.parse(bed_iter):
        for i,e2 in enumerate(g.exons):
            start2, end2 = e2
            if (start2, end2) not in exons[g.chrom]:
                exons[g.chrom].add((start2, end2))    
                exon_trees[g.chrom].insert_interval(Interval(start2, end2))
            if i > 0:
                start1, end1 = g.exons[i-1]
                if (end1, start2) not in introns[g.chrom]:
                    introns[g.chrom].add((end1, start2))    
                    intron_trees[g.chrom].insert_interval(Interval(end1, start2))
    return exon_trees, intron_trees

def get_overlap_fraction(tx, trees):
    # Determine percent overlap with annotated genes
    overlap_fraction = 0.0
    if tx.chrom in trees:
        tree = trees[tx.chrom]            
        overlap_bases = 0
        tx_length = 0    
        for e in tx.exons:  
            astart = e.start
            aend = e.end              
            bintervals = tree.find(astart, aend)
            max_exon_overlap_bases = 0
            for interval in bintervals:            
                bstart, bend = interval.start, interval.end
                overlap_start = astart if bstart < astart else bstart
                overlap_end = bend if bend < aend else aend
                max_exon_overlap_bases = max(max_exon_overlap_bases, (overlap_end - overlap_start))
            overlap_bases += max_exon_overlap_bases
            tx_length += (aend - astart)            
        overlap_fraction = overlap_bases / float(tx_length)
    return overlap_fraction

def find_transcript_cluster_stats(transcripts):
    cluster_tree = ClusterTree(0,1)     
    for id, f in enumerate(transcripts):
        cluster_tree.insert(f.start, f.end, id)
    # TODO: group by samples not lanes
    for start, end, ids in cluster_tree.getregions():
        strand_transcripts = [[], [], []]
        for id in ids:            
            strand_transcripts[transcripts[id].strand].append(transcripts[id])
        strand_densities = [0, 0, 0] 
        strand_recur = [0, 0, 0]
        for strand,txs in enumerate(strand_transcripts):
            # TODO: t.label is not the right attribute to use to group samples
            recur = len(set(t.label for t in txs))
            length = sum(t.length for t in txs) if len(txs) > 0 else 1
            density = sum(t.score for t in txs) / float(length)
            strand_recur[strand] = recur
            strand_densities[strand] = density
        total_density = (strand_densities[POS_STRAND] + 
                         strand_densities[NEG_STRAND])    
        if total_density == 0:
            # if there is no "stranded" coverage at this node, 
            # then assign coverage to the positive strand by convention
            pos_frac = 1.0
        else:
            # proportionally assign unstranded coverage based on amount of
            # plus and minus strand coverage
            pos_frac = strand_densities[POS_STRAND] / float(total_density)
        # allocate recurrence/density from unstranded transcripts
        strand_fracs = (pos_frac, 1.0 - pos_frac)
        strand_num_txs = [len(strand_transcripts[s]) for s in (POS_STRAND, NEG_STRAND, NO_STRAND)]
        mean_densities = [0, 0, 0]
        for s in (POS_STRAND, NEG_STRAND):
            strand_recur[s] += strand_fracs[s] * strand_recur[NO_STRAND]
            strand_densities[s] += strand_fracs[s] * strand_densities[NO_STRAND]
            strand_num_txs[s] += strand_fracs[s] * strand_num_txs[NO_STRAND]
            if strand_num_txs[s] == 0:
                mean_densities[s] = 0
            else:
                mean_densities[s] = strand_densities[s] / float(strand_num_txs[s])
        strand_recur[NO_STRAND] = max(strand_recur[:2])
        mean_densities[NO_STRAND] = max(mean_densities[:2])
        for id in ids:
            tx = transcripts[id]
            recur = strand_recur[tx.strand]
            mean_density = mean_densities[tx.strand]
            yield tx, recur, mean_density

def get_transcript_properties(transcripts, exon_trees, intron_trees):
    for tx_stats in find_transcript_cluster_stats(transcripts):
        tx, recur, mean_density = tx_stats
        exonic_frac = get_overlap_fraction(tx, exon_trees)
        intronic_frac = get_overlap_fraction(tx, intron_trees)
        yield tx, exonic_frac, intronic_frac, recur, mean_density 

def classify_transcripts(known_genes_file, gtf_file):
    # read known transcripts
    logging.debug("Reading known genes")
    exon_trees, intron_trees = \
        build_interval_trees(open(known_genes_file))
    # parse gtf file
    logging.debug("Reading transcripts")
    for locus_transcripts in parse_gtf(open(gtf_file)):
        for tx_stats in get_transcript_properties(locus_transcripts, exon_trees, intron_trees):
            print tx_stats

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")    
    parser = argparse.ArgumentParser()
    parser.add_argument("known_genes_file")
    parser.add_argument("gtf_file")
    options = parser.parse_args()
    classify_transcripts(options.known_genes_file, options.gtf_file)

if __name__ == '__main__': main()