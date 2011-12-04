'''
Created on Nov 3, 2011

@author: mkiyer
'''
import logging
import argparse

from assemblyline.lib.juncs import parse_raw_juncs

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument('juncs_files', nargs='+')
    args = parser.parse_args()
    juncs = set()
    for f in args.juncs_files:
        for chrom, start, end, strand in parse_raw_juncs(open(f)):
            juncs.add((chrom, start, end, strand))
    sorted_juncs = sorted(juncs)
    for chrom, start, end, strand in sorted_juncs:
        print '\t'.join(map(str, [chrom, start, end, strand]))

if __name__ == '__main__':
    main()