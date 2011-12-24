'''
Created on Dec 20, 2011

@author: mkiyer
'''
import logging
import argparse
import bisect
import collections
import pysam

from assemblyline.lib.juncs import parse_raw_juncs

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument('--dist', dest="distance", type=int, default=300)
    parser.add_argument('genome_fasta_file')
    parser.add_argument('full_junc_file')
    parser.add_argument('junc_file')
    args = parser.parse_args()
    maxdist = args.distance
    # open fasta
    fastafh = pysam.Fastafile(args.genome_fasta_file)
    # read full junc file
    splice_sites = {"+": collections.defaultdict(lambda: (set(),set())),
                    "-": collections.defaultdict(lambda: (set(),set()))}
    LEFT = 0
    RIGHT = 1
    for chrom, start, end, strand in parse_raw_juncs(open(args.full_junc_file), header=False):
        splice_sites[strand][chrom][LEFT].add(start)
        splice_sites[strand][chrom][RIGHT].add(end)
    # sort positions
    for strand, strand_sites in splice_sites.iteritems():
        for chrom, sites in strand_sites.iteritems():
            newleft = sorted(sites[LEFT])
            newright = sorted(sites[RIGHT])
            strand_sites[chrom] = (newleft, newright)
    # extract junc sequences
    f = open(args.junc_file)
    print f.next()
    for line in f:
        fields = line.strip().split('\t')
        chrom = fields[0]
        start = int(fields[1])
        end = int(fields[2])
        strand = fields[3]
        # locate nearest RIGHT before LEFT
        left_ind = bisect.bisect(splice_sites[strand][chrom][RIGHT], start)
        if left_ind == 0:
            left_pos = 0
        else:
            left_pos = splice_sites[strand][chrom][RIGHT][left_ind - 1] 
        # locate nearest LEFT after RIGHT
        right_ind = bisect.bisect_left(splice_sites[strand][chrom][LEFT], end)
        if right_ind >= len(splice_sites[strand][chrom][LEFT]):
            right_pos = 1000000000
        else:
            right_pos = splice_sites[strand][chrom][LEFT][right_ind]
        #print "START", start, "LEFT_POS", left_pos    
        #print "END", end, "RIGHT_POS", right_pos    
        # extract sequence
        left_pos = max(0, left_pos, start - maxdist)
        right_pos = min(right_pos, end + maxdist)
        left_seq = fastafh.fetch(chrom, left_pos, start)
        right_seq = fastafh.fetch(chrom, end, right_pos)
        fields.extend((left_seq, right_seq))
        print '\t'.join(map(str, fields))
        
    fastafh.close()

if __name__ == '__main__':
    main()
