'''
Created on Jan 12, 2013

@author: mkiyer
'''
import logging
import argparse
import os
import sys
import subprocess

import pysam
from assemblyline.rnaseq.lib.base import check_executable

def main():
    # setup logging
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    # command line parsing
    parser = argparse.ArgumentParser()
    parser.add_argument('--fusion-tag', dest='fusion_tag', default='XF')
    parser.add_argument('input_file')
    parser.add_argument('output_file')
    args = parser.parse_args()
    if not check_executable("samtools"):
        parser.error("samtools binary not found")
    if not os.path.exists(args.input_file):
        parser.error("Input file not found")
    fusion_tag = args.fusion_tag
    infh = pysam.Samfile(args.input_file, 'rb')
    outfh = pysam.Samfile(args.output_file, 'wb', template=infh)
    # iterate through input file
    fusion_reads = 0
    num_reads = 0
    for r in infh:
        # get tags
        tagdict = dict(r.tags)
        # write if read has fusion tag
        if fusion_tag in tagdict:
            outfh.write(r)
            fusion_reads += 1
        num_reads += 1
    outfh.close()
    infh.close()
    logging.info("Found %d fusion reads" % (fusion_reads))
    logging.info("Found %d total reads" % (num_reads))
    logging.info("Fusion read fraction: %f" % (float(fusion_reads)/num_reads))
    logging.debug("Indexing BAM file")
    retcode = subprocess.call(["samtools", "index", args.output_file])    
    logging.info("Done")
    return retcode

if __name__ == '__main__':
    sys.exit(main())