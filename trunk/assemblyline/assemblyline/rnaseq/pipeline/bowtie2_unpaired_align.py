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
    
def bowtie2_unpaired_align(bowtie2_index,
                           result_file,
                           fastq_files,
                           extra_args,
                           num_processors,
                           tmp_dir):
    num_threads = max(1, num_processors-1)
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
    parser.add_argument("bowtie2_index")
    parser.add_argument("result_file")
    parser.add_argument("fastq_files", nargs="+")
    args = parser.parse_args()
    return bowtie2_unpaired_align(args.bowtie2_index,
                                  args.result_file,
                                  args.fastq_files,
                                  args.extra_args,
                                  args.num_processors,
                                  args.tmp_dir)

if __name__ == '__main__': 
    sys.exit(main())