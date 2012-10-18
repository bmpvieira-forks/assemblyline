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

from assemblyline.lib.sampletable import LibInfo
import assemblyline.utils
utils_dir = assemblyline.utils.__path__[0]
bam_to_sam_sr_script = os.path.join(utils_dir, "bam_to_sam_sr.py")

# test imports
import HTSeq
import pysam

def run_htseq_count(bam_file, gtf_file, output_dir, library_type, 
                    sample_id, library_id):
    # bam -> sam (converts PE reads to SR)
    args = ["python", bam_to_sam_sr_script, bam_file]
    bam2sam_p = subprocess.Popen(args, stdout=subprocess.PIPE)
    # need to quantify as unstranded to remove bias
    stranded = "no"
    #if library_type == "fr-unstranded":
    #    stranded = "no"
    #else:
    #    stranded = "yes"    
    # htseq-count
    args = ["htseq-count", "-m", "union", "-s", stranded, "-", gtf_file]
    output_file = os.path.join(output_dir, "htseq_count.txt")
    outfh = open(output_file, "w")
    logging.debug("htseq-count args: %s" % (map(str, args)))
    retcode = subprocess.call(args, stdin=bam2sam_p.stdout, stdout=outfh)
    outfh.close()
    if retcode != 0:
        bam2sam_p.kill()
        if os.path.exists(output_file):
            os.remove(output_file)
        return 1, sample_id, library_id
    retcode = bam2sam_p.wait()
    if retcode != 0:
        return 1, sample_id, library_id
    # create job.done file
    open(os.path.join(output_dir, "job.done"), "w").close()
    return 0, sample_id, library_id

def run_htseq_count_task(args):
    return run_htseq_count(**args)

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
    libinfos = []
    for libinfo in LibInfo.from_file(args.library_table):
        if not libinfo.is_valid():
            logging.warning("\tskipping cohort=%s patient=%s sample=%s lib=%s lanes=%s" % 
                            (libinfo.cohort, libinfo.patient, libinfo.sample, 
                             libinfo.library, libinfo.lanes))
            continue
        libinfos.append(libinfo)
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
                        'library_type': libinfo.library_type,
                        'library_id': libinfo.library,
                        'sample_id': libinfo.sample}
            tasks.append(arg_dict)
    # use multiprocessing to parallelize
    logging.info("Running tasks")
    pool = Pool(processes=num_processes)
    result_iter = pool.imap_unordered(run_htseq_count_task, tasks)
    for retcode,sample_id,library_id in result_iter:        
        logging.debug("\tfinished sample %s library %s with return code %d" % (sample_id, library_id, retcode))
    pool.close()

    
if __name__ == '__main__':
    sys.exit(main())