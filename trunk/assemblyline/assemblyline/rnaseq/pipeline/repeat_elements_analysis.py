'''
Created on Dec 6, 2012

@author: mkiyer
'''
import argparse
import logging
import sys
import array

# local import
import pysam
# project imports
import assemblyline.rnaseq.lib.config as config

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("sam_file")
    parser.add_argument("output_file")
    args = parser.parse_args()
    samfh = pysam.Samfile(args.sam_file, "r")
    # init counts
    counts = array.array('L', (0 for x in xrange(len(samfh.references))))
    unmapped = 0
    for r in samfh:
        if r.is_unmapped:
            unmapped += 1
        else:
            counts[r.tid] += 1
    # print output
    outfh = open(args.output_file, "w")
    print >>outfh, '\t'.join(["#rname", "length", "count", "unmapped"])
    for i in xrange(len(samfh.references)):
        print >>outfh, '\t'.join(map(str, [samfh.references[i], samfh.lengths[i], counts[i], "0"]))
    print >>outfh, '\t'.join(["unmapped", "0", "0", str(unmapped)])
    return config.JOB_SUCCESS

if __name__ == '__main__': 
    sys.exit(main())