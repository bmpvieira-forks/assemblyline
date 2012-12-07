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
import assemblyline.rnaseq.lib.config as config
import assemblyline.rnaseq.pipeline
_pipeline_dir = assemblyline.rnaseq.pipeline.__path__[0]

def align_repeat_elements(bowtie2_index,
                          result_file,
                          fastq_files,
                          num_processors):
    num_threads = max(1, num_processors-1)
    # start bowtie2
    args = ["bowtie2", "-p", num_threads,
            "--local", "--very-sensitive-local", 
            "-x", bowtie2_index,
            "-U", ','.join(fastq_files)]
    logging.debug("bowtie args: %s" % (' '.join(map(str, args))))
    aln_p = subprocess.Popen(map(str,args), stdout=subprocess.PIPE)
    # process alignments
    args = [sys.executable, 
            os.path.join(_pipeline_dir, "repeat_elements_analysis.py"),
            "-", result_file]
    logging.debug("analysis args: %s" % (' '.join(map(str, args))))
    retcode = subprocess.call(map(str,args), stdin=aln_p.stdout)    
    aln_p.terminate()
    if retcode != 0:
        logging.error("repeat analysis failed")
        # cleanup output file
        if os.path.exists(result_file):
            os.remove(result_file)
        return config.JOB_ERROR
    return config.JOB_SUCCESS

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", dest="num_processors", type=int, default=1)
    parser.add_argument("bowtie2_index")
    parser.add_argument("result_file")
    parser.add_argument("fastq_files", nargs="+")
    args = parser.parse_args()
    return align_repeat_elements(args.bowtie2_index,
                                 args.result_file,
                                 args.fastq_files,
                                 args.num_processors)

if __name__ == '__main__': 
    sys.exit(main())