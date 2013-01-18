'''
Created on Dec 5, 2012

@author: mkiyer
'''
import argparse
import logging
import sys
import os
import subprocess

from assemblyline.lib.base import check_executable
import assemblyline.rnaseq.pipeline
_pipeline_dir = assemblyline.rnaseq.pipeline.__path__[0]

def bowtie2_paired_align(bowtie2_index,
                         bam_file,
                         read1_files,
                         read2_files,
                         unpaired_files,
                         extra_args,
                         num_processors,
                         tmp_dir):
    num_threads = max(1, num_processors-1)
    # align with bowtie2
    bam_prefix = os.path.splitext(bam_file)[0]
    unsorted_bam_file = os.path.join(tmp_dir, os.path.basename(bam_prefix) + ".unsorted.bam")
    args = ["bowtie2"]
    for arg in extra_args:        
        args.extend(arg.split())    
    args.extend(["-p", num_threads,
                 # TODO: problem in bowtie2 version 2.0.5 
                 # with the --no-unal flag
                 #"--no-unal",
                 "-x", bowtie2_index,
                 "-1", ','.join(read1_files),
                 "-2", ','.join(read2_files),
                 "-U", ','.join(unpaired_files)])
    args = map(str, args)
    logging.debug("bowtie2 args: %s" % (' '.join(args)))
    aln_p = subprocess.Popen(args, stdout=subprocess.PIPE)
    # remove unmapped reads from SAM
    # TODO: this will not be necessary after bug fix in bowtie2 --no-unal option
    args = ["python", 
            os.path.join(_pipeline_dir, "sam_no_unal.py"),
            "-", "-"]
    logging.debug("sam_no_unal args: %s" % str(args))
    nounal_p = subprocess.Popen(map(str,args), stdin=aln_p.stdout, 
                                stdout=subprocess.PIPE)  
    # convert sam to bam
    args = ["samtools", "view", "-bS", "-"]
    f = open(unsorted_bam_file, "wb")
    retcode1 = subprocess.call(args, stdin=nounal_p.stdout, stdout=f)     
    f.close()
    retcode2 = nounal_p.wait()
    retcode3 = aln_p.wait()
    retcode = retcode1 + retcode2 + retcode3
    if retcode != 0:
        logging.error("bowtie2 alignment failed")
        if os.path.exists(unsorted_bam_file):
            os.remove(unsorted_bam_file)
        return retcode
    # sort bam
    args = ["samtools", "sort", "-m", "2000000000",
            unsorted_bam_file, bam_prefix]
    logging.debug("samtools sort args: %s" % args)
    retcode = subprocess.call(args)
    if retcode != 0:
        logging.error("samtools sort failed")
        if os.path.exists(bam_file):
            os.remove(bam_file)
        return retcode
    bam_index = bam_file + ".bai"
    args = ["samtools", "index", bam_file]
    logging.debug("samtools index args: %s" % args)
    retcode = subprocess.call(args)
    if retcode != 0:
        logging.error("samtools index failed")
        if os.path.exists(bam_index):
            os.remove(bam_index)
        return retcode
    return 0

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", dest="num_processors", type=int, default=1)
    parser.add_argument("--tmp-dir", dest="tmp_dir", default="/tmp")
    parser.add_argument("--extra-arg", dest="extra_args", action="append", default=[])
    parser.add_argument("-1", dest="read1_files", default=None)
    parser.add_argument("-2", dest="read2_files", default=None)
    parser.add_argument("-U", dest="unpaired_files", default=None)
    parser.add_argument("bowtie2_index")
    parser.add_argument("bam_file")
    args = parser.parse_args()
    if not check_executable("bowtie2"):
        parser.error("bowtie2 binary not found")
    if not check_executable("samtools"):
        parser.error("samtools binary not found")
    read1_files = []
    read2_files = []
    unpaired_files = []
    if args.read1_files is not None:
        read1_files = args.read1_files.split(',')
        for f in read1_files:
            if not os.path.exists(f):
                parser.error("fastq file %s not found" % (f))
    if args.read2_files is not None:
        read2_files = args.read2_files.split(',')
        for f in read2_files:
            if not os.path.exists(f):
                parser.error("fastq file %s not found" % (f))
    if args.unpaired_files is not None:
        unpaired_files = args.unpaired_files.split(',')
        for f in unpaired_files:
            if not os.path.exists(f):
                parser.error("fastq file %s not found" % (f))
    if len(read1_files) != len(read2_files):
        parser.error("number of read1 files != number of read2 files")
    if not os.path.exists(args.tmp_dir):
        parser.error("tmp dir %s not found" % (args.tmp_dir))
    suffix = os.path.splitext(args.bam_file)[-1]
    if suffix != ".bam":
        parser.error("bam file %s must have .bam suffix" % (args.bam_file))    
    tmp_dir = os.path.abspath(args.tmp_dir)
    return bowtie2_paired_align(args.bowtie2_index,
                                args.bam_file,
                                read1_files,
                                read2_files,
                                unpaired_files,
                                args.extra_args,
                                args.num_processors,
                                tmp_dir)

if __name__ == '__main__': 
    sys.exit(main())