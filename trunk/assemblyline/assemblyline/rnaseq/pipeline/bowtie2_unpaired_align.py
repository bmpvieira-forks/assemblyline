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

def bowtie2_align(bowtie2_index,
                  result_file,
                  bam_file,
                  fastq_files,
                  extra_args,
                  num_processors,
                  tmp_dir):
    num_threads = max(1, num_processors-1)
    if bam_file is not None:
        # extract fastq from bam
        prefix = os.path.splitext(os.path.basename(bam_file))[0]
        prefix = os.path.join(tmp_dir, prefix)
        fastq_files = [prefix + "_1.fq", 
                       prefix + "_2.fq"]    
        args = [sys.executable,
                os.path.join(_pipeline_dir, "bam_to_fastq.py"),
                bam_file, prefix]
        logging.debug("bam_to_fastq args: %s" % (map(str, args)))
        retcode = subprocess.call(args)
        if retcode != 0:
            for f in fastq_files:
                if os.path.exists(f):
                    os.remove(f)
    if len(fastq_files) == 0:
        logging.error("No fastq files found")
        return 1
    # start bowtie2
    args = ["bowtie2"]
    for arg in extra_args:        
        args.extend(arg.split())    
    args.extend(["-p", num_threads,
                 "-x", bowtie2_index,
                 "-U", ','.join(fastq_files)])
    logging.debug("bowtie args: %s" % (' '.join(map(str, args))))
    aln_p = subprocess.Popen(map(str,args), stdout=subprocess.PIPE)
    # process alignments
    args = [sys.executable, 
            os.path.join(_pipeline_dir, "bowtie2_unpaired_count.py"),
            "-", result_file]
    logging.debug("analysis args: %s" % (' '.join(map(str, args))))
    retcode1 = subprocess.call(map(str,args), stdin=aln_p.stdout)    
    retcode2 = aln_p.wait()
    retcode = retcode1 + retcode2
    if bam_file is not None:
        for f in fastq_files:
            if os.path.exists(f):
                os.remove(f)        
    if retcode != 0:
        logging.error("bowtie2_unpaired_count.py script failed")
        # cleanup
        if os.path.exists(result_file):
            os.remove(result_file)
        return retcode    
    return 0

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", dest="num_processors", type=int, default=1)
    parser.add_argument("--tmp-dir", dest="tmp_dir", default="/tmp")
    parser.add_argument("--extra-arg", dest="extra_args", action="append", default=[])
    parser.add_argument("--bam-file", dest="bam_file", default=None)
    parser.add_argument("--fastq-file", dest="fastq_files", action="append", default=[])
    parser.add_argument("bowtie2_index")
    parser.add_argument("result_file")
    args = parser.parse_args()
    if args.bam_file is not None:
        if not os.path.exists(args.bam_file):
            parser.error("bam file %s not found" % (args.bam_file))
    return bowtie2_align(args.bowtie2_index,
                         args.result_file,
                         args.bam_file,
                         args.fastq_files,
                         args.extra_args,
                         args.num_processors,
                         args.tmp_dir)

if __name__ == '__main__': 
    sys.exit(main())