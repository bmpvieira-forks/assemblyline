'''
Created on Dec 11, 2011

@author: mkiyer
'''
import sys
import os
import logging
import argparse
import subprocess
from multiprocessing import Pool

from assemblyline.lib.sampletable import LibInfo

def run_cufflinks(bam_file,
                  output_dir,
                  gtf_file,
                  genome_fasta_file,
                  library_type,                  
                  output_label,
                  num_processors=1,
                  frag_len_mean=None,
                  frag_len_std_dev=None,
                  cufflinks_bin="cufflinks"):
    # save current dir and change to working dir
    curdir = os.getcwd()
    os.chdir(output_dir)
    # setup command line args
    args = [cufflinks_bin,
            "--no-update-check",
            "-v",
            "-p", num_processors,
            "-G", gtf_file,
            "-b", genome_fasta_file,
            "-u",
            "-L", output_label,
            "--library-type", library_type,
            "--compatible-hits-norm"]
    if frag_len_mean is not None:
        args.extend(["-m", frag_len_mean,
                     "-s", frag_len_std_dev])
    args.append(bam_file)
    log_file = "cufflinks.log"
    log_fh = open(log_file, "w")
    retcode = subprocess.call(map(str, args), stdout=log_fh, stderr=log_fh)
    log_fh.close()
    if retcode == 0:
        open("job.done", "w").close()
    os.chdir(curdir)
    return retcode, output_label

def run_cufflinks_task(args):
    return run_cufflinks(**args)

def main():
    # Command line parsing
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--output-dir', dest="output_dir", default=".")
    parser.add_argument('-p', '--num-processors', dest="num_processors", type=int, default=1)
    parser.add_argument('--ppj', dest="processors_per_job", type=int, default=1)
    parser.add_argument('--cufflinks-bin', dest="cufflinks_bin", default="cufflinks")
    parser.add_argument('genome_fasta_file')
    parser.add_argument('gtf_file')
    parser.add_argument('library_table')
    args = parser.parse_args()
    # check args
    assert args.num_processors % args.processors_per_job == 0
    num_processes = args.num_processors / args.processors_per_job
    gtf_file = os.path.abspath(args.gtf_file)
    genome_fasta_file = os.path.abspath(args.genome_fasta_file)    
    # parse library table
    logging.info("Parsing library table")
    tasks = []
    for libinfo in LibInfo.from_file(args.library_table):
        if not libinfo.is_valid():
            logging.warning("\tskipping cohort=%s patient=%s sample=%s lib=%s lanes=%s" % 
                            (libinfo.cohort, libinfo.patient, libinfo.sample, 
                             libinfo.library, libinfo.lanes))
            continue
        lib_dir = os.path.join(args.output_dir, libinfo.sample, libinfo.library)
        if not os.path.exists(lib_dir):
            logging.debug("Creating directory %s" % (lib_dir))
            os.makedirs(lib_dir)
        # check for output
        job_done_file = os.path.join(lib_dir, "job.done")
        msg = "Adding task for sample=%s library=%s" % (libinfo.sample, libinfo.library)
        if os.path.exists(job_done_file):
            logging.info("[SKIPPED] %s" % msg)
        else:
            logging.info(msg)        
            arg_dict = {'bam_file': libinfo.bam_file,
                        'output_dir': lib_dir,
                        'gtf_file': gtf_file,
                        'genome_fasta_file': genome_fasta_file,
                        'library_type': libinfo.library_type,
                        'output_label': libinfo.library,
                        'num_processors': args.processors_per_job,
                        'frag_len_mean': libinfo.frag_len_mean,
                        'frag_len_std_dev': libinfo.frag_len_std_dev,
                        'cufflinks_bin': args.cufflinks_bin} 
            tasks.append(arg_dict)
    # use multiprocessing to parallelize classification
    logging.info("Running tasks")
    pool = Pool(processes=num_processes)
    result_iter = pool.imap_unordered(run_cufflinks_task, tasks)
    for retcode,outprefix in result_iter:
        if retcode == 0:
            logging.debug("\tfinished %s" % (outprefix))
        else:
            logging.debug("\terror %s" % (outprefix))


if __name__ == '__main__':
    sys.exit(main())