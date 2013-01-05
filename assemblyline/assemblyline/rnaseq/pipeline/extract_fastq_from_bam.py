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

from assemblyline.rnaseq.lib.base import DNA_reverse_complement
import assemblyline.rnaseq.pipeline
_pipeline_dir = assemblyline.rnaseq.pipeline.__path__[0]

def parse_bam(bamfh, readnum_in_qname, remove_suffix):
    '''
    reads must be sorted by queryname 
    '''
    pe_reads = ([], [])
    num_reads = 0
    prev_qname = None    
    for r in bamfh:
        suffix = r.qname[-2:]
        # determine whether this is read1 or read2
        if readnum_in_qname:
            # use "/1" or "/2" suffix to determine read number
            if (suffix == "/1"):
                readnum = 0
            elif (suffix == "/2"):
                readnum = 1
            else:
                raise Exception("suffix /1 or /2 not found in qname")
        else:
            readnum = 1 if r.is_read2 else 0
        # optionally remove the /1 or /2 suffix
        if (remove_suffix and 
            (suffix == "/1" or suffix == "/2")):
            r.qname = r.qname[:-2]
        # if query name changes we have completely finished
        # the fragment and can reset the read data
        if (num_reads > 0) and (r.qname != prev_qname):
            yield pe_reads
            # reset state variables
            pe_reads = ([], [])
            num_reads = 0
        pe_reads[readnum].append(r)
        prev_qname = r.qname
        num_reads += 1
    if num_reads > 0:
        yield pe_reads

def remove_multihits(pe_reads):
    r1,r2 = None,None
    if len(pe_reads[0]) > 0:
        r1 = pe_reads[0][0]
    if len(pe_reads[1]) > 0:
        r2 = pe_reads[1][0]
    return r1,r2

def to_fastq(r, readnum):
    seq = DNA_reverse_complement(r.seq) if r.is_reverse else r.seq
    qual = r.qual[::-1] if r.is_reverse else r.qual 
    return "@%s/%d\n%s\n+\n%s" % (r.qname, readnum+1, seq, qual)

def bam_to_fastq(bam_file, fastq_prefix, assume_sorted, 
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
    for pe_reads in parse_bam(bamfh, readnum_in_qname, 
                              remove_suffix=True):
        r1, r2 = remove_multihits(pe_reads)
        if (r1 is not None) and (r2 is not None):
            # pair
            print >>paired_fhs[0], to_fastq(r1, 0)
            print >>paired_fhs[1], to_fastq(r2, 1)
            paired_counts += 1
        elif (r1 is not None):
            # unpaired read 1
            print >>unpaired_fhs[0], to_fastq(r1, 0)
            unpaired_counts[0] += 1
        elif (r2 is not None):
            # unpaired read 2
            print >>unpaired_fhs[1], to_fastq(r2, 1)
            unpaired_counts[1] += 1
        else:
            raise Exception("Error")
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
                        args.readnum_in_qname,
                        tmp_dir)

if __name__ == '__main__':
    sys.exit(main())
