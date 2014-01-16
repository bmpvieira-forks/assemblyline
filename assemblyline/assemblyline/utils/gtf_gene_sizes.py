'''
Created on Jan 6, 2014

@author: mkiyer
'''
import collections
import logging
import argparse

from assemblyline.lib.transcript import parse_gtf
from assemblyline.lib.bx.cluster import ClusterTree

def cluster_isoforms(transcripts):
    # cluster exons
    cluster_tree = ClusterTree(0,1)
    for t in transcripts:
        for e in t.exons:
            cluster_tree.insert(e.start, e.end, 1)
    exons = []
    for start,end,indexes in cluster_tree.getregions():
        exons.append((start,end))
    del cluster_tree
    return exons

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("gtf_file")
    args = parser.parse_args()    
    for locus_transcripts in parse_gtf(open(args.gtf_file)):
        # group transcripts by gene id
        gene_transcript_map = collections.defaultdict(lambda: [])
        for t in locus_transcripts:
            gene_transcript_map[t.attrs['gene_id']].append(t)
        # categorize genes
        for gene_id, gene_transcripts in gene_transcript_map.iteritems():
            gene_exons = cluster_isoforms(gene_transcripts)
            length = sum((e[1]-e[0]) for e in gene_exons)
            print '\t'.join([gene_id, str(length)])

if __name__ == '__main__':
    main()
