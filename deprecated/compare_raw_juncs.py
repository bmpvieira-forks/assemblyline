'''
Created on Dec 20, 2011

@author: mkiyer
'''
import logging
import argparse
import pysam

from assemblyline.lib.juncs import parse_raw_juncs

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument('junc_subtract_file')
    parser.add_argument('junc_file')
    args = parser.parse_args()
    subtract_juncs = set()
    for chrom, start, end, strand in parse_raw_juncs(open(args.junc_subtract_file), header=False):
        subtract_juncs.add((chrom, start, end, strand))
    f = open(args.junc_file)
    f.next()
    for line in f:
        fields = line.strip().split('\t')
        chrom = fields[0]
        start = int(fields[1])
        end = int(fields[2])
        strand = fields[3]
        bubba = (chrom, start, end, strand) in subtract_juncs
        fields.append(int(bubba))
        print '\t'.join(map(str, fields))

if __name__ == '__main__':
    main()
