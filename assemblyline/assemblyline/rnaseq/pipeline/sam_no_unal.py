'''
Created on Jan 17, 2013

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
    # open input file
    if args.input_file == "-":
        inmode = "r"
    else:
        if not os.path.exists(args.input_file):
            parser.error("input file not found")
        ext = os.path.splitext(args.input_file)[-1]
        if ext == ".bam":
            inmode = "rb"
        else:
            inmode = "r"
    infh = pysam.Samfile(args.input_file, inmode)
    # open output file
    if args.output_file == "-":
        outmode = "wh"
    else:
        ext = os.path.splitext(args.output_file)[-1]
        if ext == ".bam":
            outmode = "wb"
        else:
            outmode = "wh"
    outfh = pysam.Samfile(args.output_file, outmode, template=infh)
    num_reads = 0
    unmapped = 0
    for r in infh:
        if r.is_unmapped:
            unmapped += 1
            continue
        outfh.write(r)
        num_reads += 1
    logging.debug("sam_no_unal: %d unmapped reads" % (unmapped))
    logging.debug("sam_no_unal: %d total reads" % (num_reads))
    outfh.close()
    infh.close()
    return 0

if __name__ == '__main__':
    sys.exit(main())
