'''
Created on Dec 11, 2011

@author: mkiyer
'''
import sys
import os
import logging
import argparse
from multiprocessing import Pool

from assemblyline.lib.sampletable import LibInfo
from assemblyline.lib.quantify import quantify_gene_abundance

def run_quantify(bam_file, gtf_file, output_dir, sample_id, library_id):
    output_file = os.path.join(output_dir, "gene_abundance.txt")
    # ensure BAM index can be found by pysam module    
    bai_file = bam_file + ".bai"
    split_bai_file = os.path.splitext(bam_file)[0] + ".bai"
    made_symlink = False
    if not os.path.exists(bai_file):
        if not os.path.exists(split_bai_file):
            logging.error("BAM index not found")
            return 1, sample_id, library_id
        os.symlink(split_bai_file, bai_file)
        made_symlink = True
    # run abundance estimation
    try:
        quantify_gene_abundance(gtf_file, bam_file, output_file)
    except:
        logging.error("Exception caught")
        return 2, sample_id, library_id
    if made_symlink:
        #os.remove(bai_file)
        pass
    # create job.done file
    open(os.path.join(output_dir, "job.done"), "w").close()
    return 0, sample_id, library_id

def run_quantify_task(args):
    return run_quantify(**args)

def main():
    # Command line parsing
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--output-dir', dest="output_dir", default=".")
    parser.add_argument('-p', '--num-processors', dest="num_processors", type=int, default=1)
    parser.add_argument('gtf_file')
    parser.add_argument('library_table')
    args = parser.parse_args()
    # check args
    num_processes = max(1, args.num_processors)
    gtf_file = os.path.abspath(args.gtf_file)
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
                        'gtf_file': gtf_file,
                        'output_dir': lib_dir,
                        'library_id': libinfo.library,
                        'sample_id': libinfo.sample}
            tasks.append(arg_dict)
    # use multiprocessing to parallelize classification
    logging.info("Running tasks")
    pool = Pool(processes=num_processes)
    result_iter = pool.imap_unordered(run_quantify_task, tasks)
    for retcode,sample_id,library_id in result_iter:        
        logging.debug("\tfinished sample %s library %s with return code %d" % (sample_id, library_id, retcode))


if __name__ == '__main__':
    sys.exit(main())