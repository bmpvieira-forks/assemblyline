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

import assemblyline.rnaseq.pipeline
_pipeline_dir = assemblyline.rnaseq.pipeline.__path__[0]

def run_htseq_count(gtf_file, bam_file, output_file, run_as_pe, 
                    is_sorted, extra_args, tmp_dir):
    # to use paired-end alignments the bam file needs to be sorted by 
    # queryname before running htseq-count
    do_sort = run_as_pe and (not is_sorted)
    if do_sort:
        # sort bam by read name
        prefix = os.path.join(tmp_dir, os.path.splitext(os.path.basename(bam_file))[0])
        sorted_bam_prefix = prefix + ".qnamesorted"
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
        input_file = sorted_bam_file
    else:
        input_file = bam_file
    # now stream the bam file as sam for htseq-count
    if run_as_pe:
        args = ["samtools", "view", input_file]
        logging.debug("samtools view args: %s" % (map(str, args)))
        sam_p = subprocess.Popen(args, stdout=subprocess.PIPE)
    else:
        # if running as single-read need to eliminate all paired-end flags from
        # the reads 
        args = ["python", 
                os.path.join(_pipeline_dir, "bam_pe_to_sr.py"),
                input_file, "-"]
        logging.debug("bam pe to sr args: %s" % (map(str, args)))
        sam_p = subprocess.Popen(args, stdout=subprocess.PIPE)
    # run htseq-count
    args = ["htseq-count"]
    for arg in extra_args:        
        args.extend(arg.split())
    args.extend(["-", gtf_file])
    logging.debug("htseq-count args: %s" % (map(str, args)))
    outfh = open(output_file, "w")
    retcode1 = subprocess.call(args, stdin=sam_p.stdout, stdout=outfh)
    outfh.close()
    retcode2 = sam_p.wait()
    retcode = retcode1 + retcode2
    # clean up
    if do_sort and os.path.exists(sorted_bam_file):
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
    parser.add_argument("--tmp-dir", dest='tmp_dir', default="/tmp")
    parser.add_argument("--arg", dest="extra_args", action="append", default=[])
    parser.add_argument("--pe", dest='pe', action="store_true", default=False)
    parser.add_argument("--assume-sorted", dest='assume_sorted', action="store_true", default=False)
    parser.add_argument('gtf_file')
    parser.add_argument('bam_file')
    parser.add_argument('output_file')
    args = parser.parse_args()
    # check args
    if not os.path.exists(args.gtf_file):
        logging.error("gtf file %s not found" % (args.gtf_file))
        return 1
    if not os.path.exists(args.bam_file):
        logging.error("bam file %s not found" % (args.bam_file))
        return 1
    return run_htseq_count(args.gtf_file, args.bam_file, args.output_file, 
                           run_as_pe=args.pe,
                           is_sorted=args.assume_sorted,
                           extra_args=args.extra_args,
                           tmp_dir=args.tmp_dir)

if __name__ == '__main__':
    sys.exit(main())
