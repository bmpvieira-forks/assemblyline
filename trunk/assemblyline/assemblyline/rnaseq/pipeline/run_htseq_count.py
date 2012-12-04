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

import HTSeq

import assemblyline.rnaseq.pipeline
_pipeline_dir = assemblyline.rnaseq.pipeline.__path__[0]

def run_htseq_count(bam_file, gtf_file, output_file, tmp_dir, stranded):
    # sort bam by read name
    prefix = os.path.join(tmp_dir, os.path.splitext(os.path.basename(bam_file))[0])
    sorted_bam_prefix = prefix + ".srt"
    sorted_bam_file = sorted_bam_prefix + ".bam"
    args = ["samtools", "sort", "-m", "2000000000",
            "-n", bam_file, sorted_bam_prefix]
    logging.debug("samtools sort args: %s" % (map(str, args)))
    retcode = subprocess.call(args)
    if retcode != 0:
        if os.path.exists(sorted_bam_file):
            os.remove(sorted_bam_file)
        return 1
    logging.debug("sorted bam file: %s" % (bam_file))
    # convert to sam for htseq-count
    #args = [sys.executable, 
    #        os.path.join(_pipeline_dir, "fix_pe_qname.py"),
    #        sorted_bam_file, "-"]
    args = ["samtools", "view", sorted_bam_file]
    sam_p = subprocess.Popen(args, stdout=subprocess.PIPE)
    # run htseq-count
    args = ["htseq-count", "-m", "union", "-s", stranded, "-", gtf_file]
    outfh = open(output_file, "w")
    logging.debug("htseq-count args: %s" % (map(str, args)))
    retcode1 = subprocess.call(args, stdin=sam_p.stdout, stdout=outfh)
    retcode2 = sam_p.wait()
    retcode = retcode1 + retcode2
    logging.info("finished htseq-count bam file: %s" % (bam_file))
    # clean up
    outfh.close()
    if os.path.exists(sorted_bam_file):
        os.remove(sorted_bam_file)
    if retcode != 0:
        if os.path.exists(output_file):
            os.remove(output_file)
        return 1
    return 0

def main():
    # setup logging
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    # command line parsing
    parser = argparse.ArgumentParser()
    parser.add_argument('--stranded', dest="stranded", default="no")
    parser.add_argument('gtf_file')
    parser.add_argument('bam_file')
    parser.add_argument('output_file')
    parser.add_argument('tmp_dir')
    args = parser.parse_args()
    # check args
    if not os.path.exists(args.gtf_file):
        logging.error("gtf file %s not found" % (args.gtf_file))
        return 1
    if not os.path.exists(args.bam_file):
        logging.error("bam file %s not found" % (args.bam_file))
        return 1
    return run_htseq_count(args.bam_file, args.gtf_file, 
                           args.output_file, args.tmp_dir, 
                           args.stranded)

if __name__ == '__main__':
    sys.exit(main())
