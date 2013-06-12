'''
Created on Jan 22, 2013

@author: mkiyer
'''
import logging
import argparse
import sys
import os

from assemblyline.lib.transcript import parse_gtf

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("gtf_file")
    parser.add_argument("region")
    args = parser.parse_args()
    if not os.path.exists(args.gtf_file):
        parser.error("GTF file %s not found" % (args.gtf_file))
    region_chrom, startend = args.region.split(":")
    start, end = startend.split("-")
    region_start = int(start)
    region_end = int(end)
    for transcripts in parse_gtf(open(args.gtf_file)):
        for t in transcripts:
            if ((t.chrom == region_chrom) and
                (t.start < region_end) and
                (t.end > region_start)):
                features = t.to_gtf_features()
                for f in features:
                    print str(f)
    logging.debug("Done")
    return 0

if __name__ == '__main__':
    sys.exit(main())