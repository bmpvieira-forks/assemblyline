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

from assemblyline.rnaseq.lib.base import parse_reads_by_qname

def bam_rm_multihits(bam_file, output_sam_file):
    bamfh = pysam.Samfile(bam_file, "rb")
    outfh = pysam.Samfile(output_sam_file, "wh", template=bamfh)
    for reads in parse_reads_by_qname(bamfh):
        r = reads[0]
        r.is_secondary = False
        outfh.write(r)
    outfh.close()
    bamfh.close()

def bam_to_fastq(bam_file, fastq_prefix, picard_dir, tmp_dir):
    # sort bam by read name
    prefix = os.path.join(tmp_dir, os.path.splitext(os.path.basename(bam_file))[0])
    sorted_bam_prefix = prefix + ".srt"
    sorted_bam_file = sorted_bam_prefix + ".bam"
    args = ["samtools", "sort", "-n", bam_file, sorted_bam_prefix]
    logging.debug("samtools sort args: %s" % (map(str, args)))
    retcode = subprocess.call(args)
    if retcode != 0:
        logging.debug("error sorting bam file" % (bam_file))
        if os.path.exists(sorted_bam_file):
            os.remove(sorted_bam_file)
        return 1
    logging.debug("sorted bam file: %s" % (bam_file))
    # remove multihit alignments
    sorted_nomulti_sam_file = sorted_bam_prefix + ".nomulti.sam"
    logging.debug("parsing and removing multihits")
    try:
        bam_rm_multihits(sorted_bam_file, sorted_nomulti_sam_file)
    except:
        logging.error("error removing multihits")
        if os.path.exists(sorted_bam_file):
            os.remove(sorted_bam_file)
        if os.path.exists(sorted_nomulti_sam_file):
            os.remove(sorted_nomulti_sam_file)
        return 1
    # extract fastq from bam
    read1_fastq = fastq_prefix + "_1.fq"
    read2_fastq = fastq_prefix + "_2.fq"
    args = ["java", "-jar",             
            os.path.join(picard_dir, "SamToFastq.jar"),
            "INPUT=%s" % (sorted_nomulti_sam_file),
            "FASTQ=%s" % (read1_fastq),
            "SECOND_END_FASTQ=%s" % (read2_fastq),
            "RE_REVERSE=true",
            "INCLUDE_NON_PF_READS=true",
            "INCLUDE_NON_PRIMARY_ALIGNMENTS=true",
            "TMP_DIR=%s" % (tmp_dir),
            "VALIDATION_STRINGENCY=SILENT"]
    logging.debug("SamToFastq args: %s" % (map(str,args)))
    retcode = subprocess.call(args)
    if os.path.exists(sorted_bam_file):
        os.remove(sorted_bam_file)
    if os.path.exists(sorted_nomulti_sam_file):
        os.remove(sorted_nomulti_sam_file)
    if retcode != 0:
        logging.error("error while running SamToFastq.jar")
        if os.path.exists(read1_fastq):
            os.remove(read1_fastq)
        if os.path.exists(read2_fastq):
            os.remove(read2_fastq)
        return 1
    return 0

def main():
    # setup logging
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    # command line parsing
    parser = argparse.ArgumentParser()
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
                        args.picard_dir, args.tmp_dir) 

if __name__ == '__main__':
    sys.exit(main())
