'''
Created on Nov 25, 2013

@author: mkiyer
'''
import argparse
import logging
import subprocess
import collections

from assemblyline.lib.bx.cluster import ClusterTree
from assemblyline.lib.transcript import parse_gtf

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument('gtf_file')
    parser.add_argument('chrom_sizes')
    parser.add_argument("output_prefix")
    args = parser.parse_args()
    # read one locus at a time
    locus_file = args.output_prefix + '.locus.bed'
    intergenic_file = args.output_prefix + '.intergenic.bed'
    intron_file = args.output_prefix + '.intron.bed'
    locus_fileh = open(locus_file, 'w')
    introns = set()
    logging.info('Parsing transcripts by locus')
    for locus_transcripts in parse_gtf(open(args.gtf_file)):
        # find borders of locus
        locus_chrom = locus_transcripts[0].chrom
        locus_start = min(t.start for t in locus_transcripts)
        locus_end = max(t.end for t in locus_transcripts)
        print >>locus_fileh, '\t'.join([locus_chrom, str(locus_start), str(locus_end)])
        # cluster locus exons
        cluster_tree = ClusterTree(0,1)
        for t in locus_transcripts:
            # update locus         
            for e in t.exons:
                cluster_tree.insert(e.start, e.end, 1)
        exon_clusters = []
        for start,end,indexes in cluster_tree.getregions():
            exon_clusters.append((start,end))
        # get intronic regions
        e1 = exon_clusters[0]
        for j in xrange(1, len(exon_clusters)):
            e2 = exon_clusters[1]
            introns.add((locus_chrom, e1[1], e2[0]))
            e1 = e2
    locus_fileh.close()
    # write introns
    logging.info('Writing introns')
    intron_fileh = open(intron_file, 'w')
    for chrom, start,end in sorted(introns):
        print >>intron_fileh, '\t'.join([chrom, str(start), str(end)])
    intron_fileh.close()
    # take complement of locus file to find intergenic regions
    logging.info('Complementing locus intervals to find intergenic regions')
    args = ['bedtools', 'complement', 
            '-i', locus_file,
            '-g', args.chrom_sizes]
    with open(intergenic_file, 'w') as f:
        subprocess.call(args, stdout=f)
    

        
                
if __name__ == '__main__':
    main()                