'''
Created on Feb 20, 2012

@author: mkiyer

AssemblyLine: transcriptome meta-assembly from RNA-Seq

Copyright (C) 2012 Matthew Iyer

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
'''
import subprocess
import sys
import os
import logging
import argparse
from multiprocessing import Pool

import assemblyline
from assemblyline.lib.librarytable import LibraryInfo

def run_htseq_count(bam_file, gtf_file, output_dir, library_id, stranded):
    # sort bam by read name
    prefix = os.path.join(output_dir, os.path.splitext(os.path.basename(bam_file))[0])
    sorted_bam_prefix = prefix + ".srt"
    sorted_bam_file = sorted_bam_prefix + ".bam"
    args = ["samtools", "sort", "-n", bam_file, sorted_bam_prefix]
    logging.debug("samtools sort args: %s" % (map(str, args)))
    retcode = subprocess.call(args)
    if retcode != 0:
        if os.path.exists(sorted_bam_file):
            os.remove(sorted_bam_file)
        return 1, library_id
    logging.info("sorted library: %s" % (library_id))
    # convert to sam for htseq-count
    args = ["samtools", "view", sorted_bam_file]
    sam_p = subprocess.Popen(args, stdout=subprocess.PIPE)
    # run htseq-count
    args = ["htseq-count", "-m", "union", "-s", "no", "-", gtf_file]
    output_file = os.path.join(output_dir, "htseq_count.txt")
    err_file = os.path.join(output_dir, "htseq_count.stderr.txt")
    outfh = open(output_file, "w")
    errfh = open(err_file, "w")
    logging.debug("htseq-count args: %s" % (map(str, args)))
    retcode1 = subprocess.call(args, stdin=sam_p.stdout, stdout=outfh, stderr=errfh)
    retcode2 = sam_p.wait()
    retcode = retcode1 + retcode2
    logging.info("counted library: %s" % (library_id))
    # clean up
    outfh.close()
    errfh.close()
    if os.path.exists(sorted_bam_file):
        os.remove(sorted_bam_file)
    if retcode != 0:
        if os.path.exists(output_file):
            os.remove(output_file)
        return 1, library_id
    # create job.done file
    open(os.path.join(output_dir, "job.done"), "w").close()    
    return 0, library_id

def run_htseq_count_task(args):
    return run_htseq_count(**args)

def main():
    # setup logging
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    logging.info("AssemblyLine %s" % (assemblyline.__version__))
    logging.info("----------------------------------")
    # command line parsing
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--output-dir', dest="output_dir", default=".")
    parser.add_argument('-p', '--num-processors', dest="num_processors", type=int, default=1)
    parser.add_argument('--stranded', dest="stranded", action="store_true", default=False) 
    parser.add_argument('gtf_file')
    parser.add_argument('library_table')
    args = parser.parse_args()
    # check args
    num_processes = max(1, args.num_processors)
    gtf_file = os.path.abspath(args.gtf_file)
    # parse library table
    logging.info("Parsing library table")
    tasks = []
    libinfos = []
    for lib in LibraryInfo.from_file(args.library_table):
        if not lib.is_valid():
            logging.warning("\tskipping lib=%s" % (lib.library_id)) 
            continue
        libinfos.append(lib)
        lib_dir = os.path.join(args.output_dir, lib.library_id)
        if not os.path.exists(lib_dir):
            logging.debug("Creating directory %s" % (lib_dir))
            os.makedirs(lib_dir)
        # check for output
        job_done_file = os.path.join(lib_dir, "job.done")
        msg = "Adding task for library=%s" % (lib.library_id)
        if os.path.exists(job_done_file):
            logging.info("[SKIPPED] %s" % msg)
        else:
            logging.info(msg)
            arg_dict = {'bam_file': lib.bam_file,
                        'gtf_file': gtf_file,
                        'output_dir': lib_dir,
                        'library_id': lib.library_id,
                        'stranded': args.stranded}
            tasks.append(arg_dict)
    # use multiprocessing to parallelize
    logging.info("Running tasks")
    pool = Pool(processes=num_processes)
    result_iter = pool.imap_unordered(run_htseq_count_task, tasks)
    for retcode,library_id in result_iter:        
        logging.debug("\tfinished library %s with return code %d" % (library_id, retcode))
    pool.close()
    pool.join()
    logging.info("Done")

if __name__ == '__main__':
    sys.exit(main())