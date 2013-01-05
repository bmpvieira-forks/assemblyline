'''
Created on Dec 4, 2012

@author: mkiyer
'''
import argparse
import logging
import os
import sys

import pysam

from assemblyline.rnaseq.lib.base import DNA_reverse_complement, make_pe_files

def bam_to_fastq(input_file, fastq_prefix):
    fastq_files = make_pe_files(fastq_prefix,'.fq')
    fastq_fhs = [open(f, "w") for f in fastq_files]
    if (input_file == "-" or os.path.splitext(input_file)[-1] == ".sam"): 
        mode = "r"
    else:
        mode = "rb"
    infh = pysam.Samfile(input_file, mode)
    for r in infh:
        seq = DNA_reverse_complement(r.seq) if r.is_reverse else r.seq
        qual = r.qual[::-1] if r.is_reverse else r.qual 
        readnum = 1 if r.is_read2 else 0
        print >>fastq_fhs[readnum], "@%s/%d\n%s\n+\n%s" % (r.qname, readnum+1, seq, qual)
    infh.close()
    for f in fastq_fhs:
        f.close()

def main():
    # setup logging
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    # command line parsing
    parser = argparse.ArgumentParser()
    parser.add_argument('input_file')
    parser.add_argument('fastq_prefix')
    args = parser.parse_args()
    return bam_to_fastq(args.input_file, args.fastq_prefix)

if __name__ == '__main__':
    sys.exit(main())