'''
Created on Dec 2, 2012

@author: mkiyer
'''
import logging
import argparse
import os
import sys
import subprocess

import pysam

from assemblyline.rnaseq.lib.libtable import FRAGMENT_LAYOUTS
import assemblyline.rnaseq.pipeline
_pipeline_dir = assemblyline.rnaseq.pipeline.__path__[0]

def bam_to_fastq(bam_file, fastq_prefix, tmp_dir,
                 fragment_layout, readnum_in_qname, do_sort):
    # sort bam by queryname
    prefix = os.path.join(tmp_dir, os.path.splitext(os.path.basename(bam_file))[0])
    sorted_bam_prefix = prefix + '.qnamesorted'
    sorted_bam_file = prefix + '.qnamesorted.bam'
    if do_sort:
        logging.info("sorting input file by queryname")
        args = ["samtools", "sort", "-m", "2000000000", 
                "-n", bam_file, sorted_bam_prefix]
        logging.debug("bam sort args: %s" % (map(str, args)))    
        retcode = subprocess.call(args)
        if retcode != 0:
            logging.error("bam sort failed")
            if os.path.exists(sorted_bam_file):
                os.remove(sorted_bam_file)
            return 1
        logging.debug("bam sort done")
    # remove multihit alignments
    logging.info("extracting individual reads")
    args = [sys.executable, 
            os.path.join(_pipeline_dir, "bam_rm_multihits.py"),
            "--fragment-layout", fragment_layout]
    if readnum_in_qname:
        args.append("--readnum-in-qname")
    args.extend([sorted_bam_file, "-"])
    logging.debug("bam_rm_multihits args: %s" % (map(str, args)))
    bam_p = subprocess.Popen(args, stdout=subprocess.PIPE)
    # extract fastq from bam
    fastq_files = [fastq_prefix + "_1.fq", 
                   fastq_prefix + "_2.fq"]    
    args = [sys.executable,
            os.path.join(_pipeline_dir, "bam_to_fastq.py"),
            "-", fastq_prefix]
    logging.debug("bam_to_fastq args: %s" % (map(str, args)))
    retcode1 = subprocess.call(args, stdin=bam_p.stdout)
    retcode2 = bam_p.wait()
    retcode = retcode1 + retcode2
    if retcode != 0:
        logging.error("error extracting reads for fastq file")
        for f in fastq_files:
            if os.path.exists(f):
                os.remove(f)
    if do_sort:
        if os.path.exists(sorted_bam_file):
            os.remove(sorted_bam_file)
    return retcode

def main():
    # setup logging
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    # command line parsing
    parser = argparse.ArgumentParser()
    parser.add_argument('--nosort', dest="do_sort", action="store_const", const=False, default=True)
    parser.add_argument('--fragment-layout', dest="fragment_layout", default="paired")
    parser.add_argument('--readnum-in-qname', dest="readnum_in_qname", action="store_true", default=False)
    parser.add_argument('--tmp-dir', dest="tmp_dir", default="/tmp")
    parser.add_argument('bam_file')
    parser.add_argument('fastq_prefix')
    args = parser.parse_args()
    if args.fragment_layout not in FRAGMENT_LAYOUTS:
        parser.error("Invalid fragment layout '%s'" % (args.fragment_layout))
    if not os.path.exists(args.bam_file):
        parser.error("Input bam file not found")
    return bam_to_fastq(args.bam_file, args.fastq_prefix, 
                        args.tmp_dir,
                        args.fragment_layout,
                        args.readnum_in_qname,
                        args.do_sort)

if __name__ == '__main__':
    sys.exit(main())
