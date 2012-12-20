'''
Created on Dec 14, 2012

@author: mkiyer
'''
import logging
import argparse
import os
import sys

import pysam

def main():
    # setup logging
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    # command line parsing
    parser = argparse.ArgumentParser()
    parser.add_argument('input_file')
    parser.add_argument('output_file')
    args = parser.parse_args()
    if not os.path.exists(args.input_file):
        parser.error("Input file not found")
    # open input file
    ext = os.path.splitext(args.input_file)[-1]
    if ext == ".sam":
        mode = "r"
    elif ext == ".bam":
        mode = "rb"
    else:
        parser.error("Unrecognized input file extension '%s' (must be .sam or .bam)" % (ext))
    infh = pysam.Samfile(args.input_file, mode)
    # open output file
    ext = os.path.splitext(args.output_file)[-1]
    if (args.output_file == "-" or ext == ".sam"):
        mode = "wh"
    elif ext == ".bam":
        mode = "wb"
    else:
        parser.error("Unrecognized output file extension '%s' (must be .sam or .bam)" % (ext))
    outfh = pysam.Samfile(args.output_file, mode, template=infh)
    num_reads = 0
    for r in infh:
        # erase 'is_paired', 'is_proper_pair', 'mate_is_unmapped', 'mate_is_reverse', 'is_read2'
        r.flag = r.flag & 0b1111111100010100
        # reset mate reference tid
        r.rnext = -1
        outfh.write(r)
        num_reads += 1
    logging.debug("bam_pe_to_sr: %d total reads" % (num_reads))
    outfh.close()
    infh.close()
    return 0

if __name__ == '__main__':
    sys.exit(main())
