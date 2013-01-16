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

from assemblyline.rnaseq.lib.base import DNA_reverse_complement, parse_sam, remove_multihits
import assemblyline.rnaseq.pipeline
_pipeline_dir = assemblyline.rnaseq.pipeline.__path__[0]

def to_fastq(r, readnum):
    seq = DNA_reverse_complement(r.seq) if r.is_reverse else r.seq
    qual = r.qual[::-1] if r.is_reverse else r.qual 
    return "@%s/%d\n%s\n+\n%s" % (r.qname, readnum+1, seq, qual)

def bam_to_fastq(bam_file, fastq_prefix, assume_sorted, has_unpaired,
                 readnum_in_qname, tmp_dir): 
    # get sort order of bam file
    sort_order = None
    f = pysam.Samfile(bam_file, 'rb')
    if 'HD' in f.header:
        hd_dict = f.header['HD']
        if 'SO' in hd_dict:
            sort_order = hd_dict['SO']
    f.close()
    # decide whether to sort bam file
    if assume_sorted or (sort_order == "queryname"):
        sorted_bam_file = bam_file
        keep_sorted_bam = True
    else:
        # sort bam by queryname
        prefix = os.path.join(tmp_dir, os.path.splitext(os.path.basename(bam_file))[0])
        sorted_bam_prefix = prefix + '.qnamesorted'
        sorted_bam_file = prefix + '.qnamesorted.bam'
        keep_sorted_bam = False
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
    # extract paired and unpaired reads
    # setup output files
    paired_files = []
    paired_fhs = []
    paired_counts = 0
    unpaired_files = []
    unpaired_fhs = []
    unpaired_counts = [0,0]
    for readnum in (1,2):
        filename = "%s.paired.%d.fq" % (fastq_prefix, readnum)
        paired_files.append(filename)
        paired_fhs.append(open(filename, 'w'))
        filename = "%s.unpaired.%d.fq" % (fastq_prefix, readnum)
        unpaired_files.append(filename)
        unpaired_fhs.append(open(filename, 'w'))
    # parse bam
    logging.debug("writing fastq files")
    bamfh = pysam.Samfile(sorted_bam_file, 'rb')
    for pe_reads in parse_sam(bamfh, readnum_in_qname, 
                              remove_suffix=True):
        r1, r2 = remove_multihits(pe_reads)
        if (r1 is None) and (r2 is None):
            raise Exception("Error parsing SAM file")
        elif (r1 is not None) and (r2 is not None):
            # pair
            print >>paired_fhs[0], to_fastq(r1, 0)
            print >>paired_fhs[1], to_fastq(r2, 1)
            paired_counts += 1
        elif has_unpaired:
            if (r1 is not None):
                # unpaired read 1
                print >>unpaired_fhs[0], to_fastq(r1, 0)
                unpaired_counts[0] += 1
            elif (r2 is not None):
                # unpaired read 2
                print >>unpaired_fhs[1], to_fastq(r2, 1)
                unpaired_counts[1] += 1
        else:
            raise Exception("Unpaired reads found when has_unpaired=False")
    bamfh.close()
    logging.debug("found %d paired reads" % (paired_counts))
    logging.debug("found %d unpaired read1" % (unpaired_counts[0]))
    logging.debug("found %d unpaired read2" % (unpaired_counts[1]))
    # clean up
    for readnum in xrange(2):
        paired_fhs[readnum].close()
        unpaired_fhs[readnum].close()
    if not keep_sorted_bam:
        if os.path.exists(sorted_bam_file):
            os.remove(sorted_bam_file)
    return 0

def main():
    # setup logging
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    # command line parsing
    parser = argparse.ArgumentParser()
    parser.add_argument('--assume-sorted', dest="assume_sorted",
                        action="store_true", default=False)
    parser.add_argument('--readnum-in-qname', dest="readnum_in_qname", 
                        action="store_true", default=False)
    parser.add_argument('--no-unpaired', dest="has_unpaired", 
                        action="store_false", default=True,
                        help="set to true if you are sure there are no "
                        "unpaired reads in the bam file. this makes error "
                        "checking more stringent")
    parser.add_argument('--tmp-dir', dest="tmp_dir", default="/tmp")
    parser.add_argument('bam_file')
    parser.add_argument('fastq_prefix')
    args = parser.parse_args()
    if not os.path.exists(args.bam_file):
        parser.error("Input bam file not found")
    if not os.path.exists(args.tmp_dir):
        parser.error("tmp dir %s not found" % (args.tmp_dir))
    tmp_dir = os.path.abspath(args.tmp_dir)
    return bam_to_fastq(args.bam_file, 
                        args.fastq_prefix, 
                        args.assume_sorted,
                        args.has_unpaired,
                        args.readnum_in_qname,
                        tmp_dir)

if __name__ == '__main__':
    sys.exit(main())
