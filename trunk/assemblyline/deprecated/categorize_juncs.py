'''
Created on Nov 22, 2011

@author: mkiyer
'''
import argparse
import logging
import os
import sys
import collections

from assemblyline.lib.juncs import parse_raw_juncs

def main():
    #
    # Command line parsing
    # 
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument('known_raw_juncs_file')
    parser.add_argument('predicted_juncs_file')
    args = parser.parse_args()
    #
    # Read known juncs file
    #
    logging.debug("Reading list of known junctions")
    known_juncs = set()
    known_left_splice_sites = set()
    known_right_splice_sites = set()    
    for chrom, start, end, strand in parse_raw_juncs(open(args.known_raw_juncs_file)):
        known_juncs.add((chrom, start, end, strand))
        known_left_splice_sites.add((chrom, start, strand))
        known_right_splice_sites.add((chrom, end, strand))        
    #
    # Read predicted juncs file
    #
    logging.debug("Categorizing predicted junctions")
    counts = collections.defaultdict(lambda: 0)
    for chrom, start, end, strand in parse_raw_juncs(open(args.predicted_juncs_file)):
        is_known = int((chrom, start, end, strand) in known_juncs)
        is_left_known = (chrom, start, strand) in known_left_splice_sites
        is_right_known = (chrom, end, strand) in known_right_splice_sites
        if is_known:
            counts["a"] += 1
        elif (is_left_known and is_right_known):
            counts["kk"] += 1
        elif (is_left_known or is_right_known):
            counts["kn"] += 1
        else:
            counts["nn"] += 1    
    logging.info("Results")
    for k,v in counts.iteritems():
        logging.info("\t%s=%d" % (k,v))
    
if __name__ == '__main__':
    sys.exit(main())
