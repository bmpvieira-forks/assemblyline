'''
Created on Dec 5, 2012

@author: mkiyer
'''
import argparse
import logging
import sys
import os
import subprocess

# project imports
import assemblyline.rnaseq.pipeline
_pipeline_dir = assemblyline.rnaseq.pipeline.__path__[0]

def bowtie2_paired_align(bowtie2_index,
                         fastq_files,
                         bam_file,
                         max_frag_size,
                         extra_args,
                         num_processors,
                         tmp_dir):
    num_threads = max(1, num_processors-1)
    bam_prefix = os.path.splitext(bam_file)[0]
    unsorted_bam_file = os.path.join(tmp_dir, os.path.basename(bam_prefix) + ".bam")
    # start bowtie2
    args = ["bowtie2"]
    for arg in extra_args:        
        args.extend(arg.split())    
    args.extend(["-p", num_threads,
                 "-I", 0,
                 "-X", max_frag_size,
                 "--no-unal",
                 "--no-mixed",
                 "--no-discordant",
                 "-x", bowtie2_index,
                 "-1", fastq_files[0],
                 "-2", fastq_files[1]])
    args = map(str, args)
    logging.debug("bowtie2 args: %s" % (' '.join(args)))
    aln_p = subprocess.Popen(args, stdout=subprocess.PIPE)
    # convert sam to bam
    args = ["samtools", "view", "-bS", "-"]
    f = open(unsorted_bam_file, "wb")
    retcode1 = subprocess.call(args, stdin=aln_p.stdout, stdout=f)     
    f.close()
    retcode2 = aln_p.wait()
    retcode = retcode1 + retcode2
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
    parser.add_argument("--max-frag-size", dest="max_frag_size", type=int, default=1000)
    parser.add_argument("bowtie2_index")
    parser.add_argument("bam_file")
    parser.add_argument("fastq_files", nargs="+")
    args = parser.parse_args()
    for filename in args.fastq_files:
        if not os.path.exists(filename):
            parser.error("fastq file %s not found" % (filename))
    bowtie2_paired_align(args.bowtie2_index,
                         args.fastq_files,
                         args.bam_file,
                         args.max_frag_size,
                         args.extra_args,
                         args.num_processors,
                         args.tmp_dir)

if __name__ == '__main__': 
    sys.exit(main())