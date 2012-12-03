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

def parse_pe_reads(bamfh, fix_qname=True):
    pe_reads = ([], [])
    # reads must be sorted by qname
    num_reads = 0
    prev_qname = None
    for read in bamfh:
        # get read attributes
        qname = read.qname
        if fix_qname:
            suffix = qname[-2:]
            if (suffix == "/1" or suffix == "/2"):
                qname = qname[:-2]
        readnum = 1 if read.is_read2 else 0
        # if query name changes we have completely finished
        # the fragment and can reset the read data
        if num_reads > 0 and qname != prev_qname:
            yield pe_reads
            # reset state variables
            pe_reads = ([], [])
            num_reads = 0
        pe_reads[readnum].append(read)
        prev_qname = qname
        num_reads += 1
    if num_reads > 0:
        yield pe_reads
        
def bam_pe_rm_multihits(bam_file, output_bam_file):
    bamfh = pysam.Samfile(bam_file, "rb")
    outfh = pysam.Samfile(output_bam_file, "wb", template=bamfh)
    num_frags = 0
    for pe_reads in parse_pe_reads(bamfh):
        n1 = len(pe_reads[0])
        n2 = len(pe_reads[1])
        if (n1 == 0 or n2 == 0):
            continue
        r1 = pe_reads[0][0]
        r2 = pe_reads[1][0]
        r1.is_secondary = False
        r2.is_secondary = False
        outfh.write(r1)
        outfh.write(r2)
        num_frags += 1
    outfh.close()
    bamfh.close()
    logging.debug("Found %d fragments" % (num_frags))

def parse_sr_reads(bamfh, fix_qname=True):
    reads = []
    for read in bamfh:        
        if len(reads) > 0 and read.qname != reads[-1].qname:
            yield reads
            reads = []
        reads.append(read)
    if len(reads) > 0:
        yield reads

def bam_sr_rm_multihits(bam_file, output_bam_file):
    bamfh = pysam.Samfile(bam_file, "rb")
    outfh = pysam.Samfile(output_bam_file, "wb", template=bamfh)
    num_frags = 0
    for reads in parse_sr_reads(bamfh):
        r = reads[0]
        r.is_secondary = False
        outfh.write(r)
        num_frags += 1
    outfh.close()
    bamfh.close()
    logging.debug("Found %d fragments" % (num_frags))

def bam_to_fastq(bam_file, fastq_prefix, picard_dir, tmp_dir,
                 fragment_layout):
    # sort bam by queryname
    prefix = os.path.join(tmp_dir, os.path.splitext(os.path.basename(bam_file))[0])
    sorted_bam_prefix = prefix + '.qnamesorted'
    sorted_bam_file = prefix + '.qnamesorted.bam'
    args = ["samtools", "sort", "-m", "2000000000", 
            "-n", bam_file, sorted_bam_prefix]
    logging.debug("bam sort args: %s" % (map(str, args)))    
#    args = ["java", "-Xmx2g", "-jar", 
#            os.path.join(picard_dir, "SortSam.jar"),
#            "INPUT=%s" % (bam_file),
#            "OUTPUT=%s" % (sorted_bam_file),
#            "SO=queryname",
#            "VALIDATION_STRINGENCY=SILENT",
#            "TMP_DIR=%s" % (tmp_dir),
#            "MAX_RECORDS_IN_RAM=1000000"]
    retcode = subprocess.call(args)
    if retcode != 0:
        logging.error("bam sort failed")
        if os.path.exists(sorted_bam_file):
            os.remove(sorted_bam_file)
        return 1
    logging.debug("bam sort done")
    # remove multihit alignments
    sorted_nomulti_bam_file = prefix + ".qnamesorted.nomulti.bam"
    logging.debug("Removing multihits")
    try:
        if fragment_layout == "paired":
            bam_pe_rm_multihits(sorted_bam_file, sorted_nomulti_bam_file)
        else:
            bam_sr_rm_multihits(sorted_bam_file, sorted_nomulti_bam_file)
    except:
        logging.error("Removing multihits failed")
        if os.path.exists(sorted_bam_file):
            os.remove(sorted_bam_file)
        if os.path.exists(sorted_nomulti_bam_file):
            os.remove(sorted_nomulti_bam_file)
        return 1
    logging.debug("Removing multihits done")    
    # extract fastq from bam
    read1_fastq = fastq_prefix + "_1.fq"
    read2_fastq = fastq_prefix + "_2.fq"
    args = ["java", "-Xmx4g", "-jar",             
            os.path.join(picard_dir, "SamToFastq.jar"),
            "RE_REVERSE=true",
            "INCLUDE_NON_PF_READS=true",
            "TMP_DIR=%s" % (tmp_dir),
            "VALIDATION_STRINGENCY=SILENT",
            "INPUT=%s" % (sorted_nomulti_bam_file),    
            "FASTQ=%s" % (read1_fastq),
            "SECOND_END_FASTQ=%s" % (read2_fastq)]
    logging.debug("SamToFastq args: %s" % (map(str, args)))
    retcode = subprocess.call(args)
    if retcode != 0:
        logging.error("error while running SamToFastq.jar")
        if os.path.exists(read1_fastq):
            os.remove(read1_fastq)
        if os.path.exists(read2_fastq):
            os.remove(read2_fastq)
    if os.path.exists(sorted_bam_file):
        os.remove(sorted_bam_file)
    if os.path.exists(sorted_nomulti_bam_file):
        os.remove(sorted_nomulti_bam_file)
    return retcode

def main():
    # setup logging
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    # command line parsing
    parser = argparse.ArgumentParser()
    parser.add_argument('--fragment-layout', dest="fragment_layout", default="paired")
    parser.add_argument('bam_file')
    parser.add_argument('fastq_prefix')
    parser.add_argument('picard_dir')
    parser.add_argument('tmp_dir')
    args = parser.parse_args()
    if not os.path.exists(args.bam_file):
        parser.error("Input bam file not found")
    if not os.path.exists(os.path.join(args.picard_dir, 'SamToFastq.jar')):
        parser.error("Picard SamToFastq.jar file not found")
    return bam_to_fastq(args.bam_file, args.fastq_prefix, 
                        args.picard_dir, args.tmp_dir,
                        args.fragment_layout)

if __name__ == '__main__':
    sys.exit(main())
