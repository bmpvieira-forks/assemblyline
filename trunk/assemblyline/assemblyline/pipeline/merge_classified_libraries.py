'''
Created on Sep 28, 2012

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
import sys
import os
import argparse
import logging

import assemblyline
import assemblyline.lib.config as config
from assemblyline.lib.gtf import merge_sort_gtf_files
from assemblyline.lib.base import CategoryCounts

def merge_transcripts(results):
    # read library category statistics
    counts_list = list(CategoryCounts.from_file(results.category_counts_file))
    library_ids = []
    expressed_gtf_files = []
    background_gtf_files = []
    for countsobj in counts_list:
        library_id = countsobj.library_id
        library_ids.append(library_id)
        prefix = os.path.join(results.classify_dir, library_id)
        expressed_gtf_files.append(prefix + ".expr.gtf")
        background_gtf_files.append(prefix + ".bkgd.gtf")
    # make a classification report
    fileh = open(results.classify_report_file, 'w')
    header_fields = ["library_id", "category", "type", "auc", "cutoff",
                     "train.tp", "train.fp", "train.fn", "train.tn",
                     "train.sens", "train.spec", "train.balacc",
                     "test.tp", "test.fp", "test.fn", "test.tn",
                     "test.sens", "test.spec", "test.balacc"]
    print >>fileh, '\t'.join(header_fields)
    for library_id in library_ids:
        prefix = os.path.join(results.classify_dir, library_id)
        intergenic_perf_file = prefix + ".intergenic.perf.txt"
        intronic_perf_file = prefix + ".intronic.perf.txt"
        input_fileh = open(intergenic_perf_file)
        input_fileh.next()
        for line in input_fileh:
            fields = [library_id, "intergenic"] + line.strip().split('\t')
            print >>fileh, '\t'.join(fields)
        input_fileh.close()
        input_fileh = open(intronic_perf_file)
        input_fileh.next()        
        for line in input_fileh:
            fields = [library_id, "intronic"] + line.strip().split('\t')
            print >>fileh, '\t'.join(fields)
        input_fileh.close()
    fileh.close()
    # add reference gtf file
    expressed_gtf_files.append(results.ref_gtf_file)
    background_gtf_files.append(results.ref_gtf_file)
    # merge sort gtf files
    logging.info("Merging and sorting expressed GTF files")
    merge_sort_gtf_files(expressed_gtf_files, 
                         results.expressed_gtf_file, 
                         tmp_dir=results.tmp_dir)
    logging.info("Merging and sorting background GTF files")
    merge_sort_gtf_files(background_gtf_files, 
                         results.background_gtf_file, 
                         tmp_dir=results.tmp_dir)
    logging.info("Done")
    return 0

def main():
    # parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("-v", "--verbose", action="store_true", 
                        dest="verbose", default=False)
    parser.add_argument("run_dir")
    args = parser.parse_args()
    # check command line parameters
    # check command line parameters
    if not os.path.exists(args.run_dir):
        parser.error("Run directory %s not found" % (args.run_dir))
    # set logging level
    if args.verbose:
        level = logging.DEBUG
    else:
        level = logging.INFO
    logging.basicConfig(level=level,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    logging.info("AssemblyLine %s" % (assemblyline.__version__))
    logging.info("----------------------------------")   
    # show parameters
    logging.info("Parameters:")
    logging.info("run directory:    %s" % (args.run_dir))
    # setup results
    results = config.AssemblylineResults(args.run_dir)
    return merge_transcripts(results) 
    logging.info("Done")
    return 0

if __name__ == '__main__':
    sys.exit(main())