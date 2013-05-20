'''
Created on Oct 17, 2012

@author: mkiyer
'''
import argparse
import logging
import os
import sys
import collections

def parse_raw_juncs(line_iter, header=False):
    """
    generator to parse lines in the tophat "raw juncs" format
    """
    if header:
        line_iter.next()
    for line in line_iter:
        if line.startswith("#"):
            continue
        fields = line.strip().split('\t')
        chrom = fields[0]
        start = int(fields[1])
        end = int(fields[2])
        strand = fields[3]
        yield (chrom, start, end, strand)

def parse_bed_juncs(line_iter):
    for line in line_iter:
        if not line:
            continue
        if line.startswith("#") or line.startswith("track"):
            continue
        fields = line.strip().split('\t')
        if len(fields) < 12:
            continue
        chrom = fields[0]
        start = int(fields[1])
        strand = fields[5]
        block_sizes = map(int, fields[10].split(','))
        block_starts = map(int, fields[11].split(','))
        for i in xrange(1,len(block_sizes)):
            if not block_sizes[i]:
                break
            left_pos = start + block_starts[i-1] + block_sizes[i-1] - 1
            right_pos = start + block_starts[i]
            yield (chrom, left_pos, right_pos, strand)

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument('ref_raw_juncs_file')    
    parser.add_argument('bed_file')
    args = parser.parse_args()     
    # parse reference juncs
    ref_juncs = collections.defaultdict(lambda: set())
    num_ref_juncs = 0
    for j in parse_raw_juncs(open(args.ref_raw_juncs_file)):
        ref_juncs[j[0]].add(j)
        num_ref_juncs += 1
    logging.info("Parsed %d reference junctions" % (num_ref_juncs))
    # parse test juncs
    known = 0
    no_chrom = 0
    novel = 0
    num_bed_juncs = 0
    for j in parse_bed_juncs(open(args.bed_file)):
        num_bed_juncs += 1
        if j[0] not in ref_juncs:
            no_chrom += 1
            continue
        chrom_juncs = ref_juncs[j[0]]
        if j in chrom_juncs:
            known += 1
        else:
            novel += 1
    print 'novel', novel
    print 'known', known
    print 'no chrom', no_chrom
    print 'total', num_bed_juncs

if __name__ == '__main__':
    sys.exit(main())
