'''
Created on Jun 9, 2012

@author: mkiyer
'''
import sys
import os
import argparse
import logging
import collections

from assemblyline.lib.bx.cluster import ClusterTree
from assemblyline.lib.transcript import parse_gtf, NO_STRAND, POS_STRAND, NEG_STRAND, strand_int_to_str

def parse_gtf_tss(gtf_file):
    tss_dict = {}
    for locus_transcripts in parse_gtf(open(gtf_file)):
        for t in locus_transcripts:            
            if t.strand == NO_STRAND:
                continue
            tss_id = t.attrs['tss_id']
            if tss_id in tss_dict:
                continue
            tss_dict[tss_id] = (t.chrom, t.strand, t.start, t.end)
    return tss_dict

def cluster_tss(tss_dict, chrom_sizes, upstream=0, downstream=0):
    trees = collections.defaultdict(lambda: ClusterTree(0,1))
    tss_ids = sorted(tss_dict)
    for i,tss_id in enumerate(tss_ids):
        chrom, strand, tstart, tend = tss_dict[tss_id]
        if strand == NEG_STRAND:
            pstart = tend - 1 - downstream
            pstart = max(tstart, pstart)
            pend = tend - 1 + upstream
            pend = min(pend, chrom_sizes[chrom])
        else:
            pstart = tstart - upstream
            pstart = max(0, pstart)
            pend = tstart + downstream
            pend = min(tend, pend)
        trees[chrom].insert(pstart, pend, i)
    # get clusters
    cur_p_id = 1
    for chrom in sorted(trees):
        tree = trees[chrom]
        for start, end, indexes in tree.getregions():
            p_id = 'P%07d' % (cur_p_id)
            for i in indexes:
                tss_id = tss_ids[i]
                yield tss_id, p_id, chrom, start, end
            cur_p_id += 1

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    # parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument('-u', type=int, dest='upstream', default=1000)
    parser.add_argument('-d', type=int, dest='downstream', default=0)
    parser.add_argument('gtf_file')
    parser.add_argument('chrom_sizes')
    args = parser.parse_args()
    upstream = args.upstream
    downstream = args.downstream
    chrom_sizes_file = args.chrom_sizes
    gtf_file = args.gtf_file
    # check command line parameters
    if not os.path.exists(args.gtf_file):
        parser.error("GTF file %s not found" % (args.gtf_file))
    logging.debug('Reading chrom size info')
    chrom_sizes = {}
    with open(chrom_sizes_file) as fileh:
        for line in fileh:
            fields = line.strip().split('\t')
            chrom_sizes[fields[0]] = int(fields[1])
    # parse
    logging.debug('Parsing TSS from GTF')
    tss_dict = parse_gtf_tss(gtf_file)
    logging.debug('Clustering TSS using upstream %dbp and downstream %dbp' % (upstream, downstream))
    print '\t'.join(['tss_id', 'p_id', 'p_chrom', 'p_start', 'p_end'])
    for fields in cluster_tss(tss_dict, chrom_sizes, upstream, downstream):
        print '\t'.join(map(str, fields))

if __name__ == '__main__':
    sys.exit(main())
