'''
Created on Nov 30, 2011

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
import logging
import argparse
import os
import subprocess
import sys
import shutil
import multiprocessing

import assemblyline
from assemblyline.lib.transcript_parser import parse_gtf
from assemblyline.lib.base import GTFAttr
from assemblyline.lib.batch_sort import batch_sort

from annotate_transcripts import LIB_COUNTS_FILE, LibCounts, CategoryInfo, \
    CATEGORIES, category_int_to_str, ANNOTATED, CATEGORY, MEAN_SCORE, \
    MEAN_RECURRENCE

# R script to call for classifying transcripts
_module_dir = assemblyline.__path__[0]
R_SCRIPT = os.path.join(_module_dir, "lib", "classify_transcripts.R")
    
def get_classification_header():
    header = ["chrom", "start", "library_id", "t_id", "annotated", 
              "category", "length", "num_exons", "score", 
              "mean_score", "mean_recurrence"]
    return header

def get_classification_fields(t):
    # setup list of annotation fields
    fields = [t.chrom,
              t.start,
              t.attrs[GTFAttr.LIBRARY_ID],
              t.attrs[GTFAttr.TRANSCRIPT_ID],
              int(t.attrs[ANNOTATED]),
              int(t.attrs[CATEGORY]),
              t.length,
              len(t.exons),
              t.score,
              t.attrs[MEAN_SCORE],
              t.attrs[MEAN_RECURRENCE]]
    return fields

def get_classification_result_header():
    header = get_classification_header()
    header.append("pred")
    return header

def run_classification_r_script(args):
    library_id, filename = args
    logging.debug("\t[STARTED]  library_id='%s'" % (library_id))
    # setup classification result files on per-library basis
    prefix = os.path.splitext(filename)[0]
    logfile = prefix + ".classify.log"
    # run R script to do classification
    logfh = open(logfile, "w")
    retcode = subprocess.call(["Rscript", R_SCRIPT, filename, prefix], 
                              stdout=logfh, stderr=logfh)
    logfh.close()
    if retcode != 0:
        logging.error("\t[FAILED]   library_id='%s'" % (library_id))
        return retcode, prefix
    # store cutoff
    logging.debug("\t[FINISHED] library_id='%s'" % (library_id))
    return retcode, prefix

def classify_category(cinfo, tasks, num_processes, tmp_dir):
    # open output file
    result_fh = open(cinfo.ctree_file, "w")
    print >>result_fh, '\t'.join(get_classification_result_header())
    # use multiprocessing to parallelize
    pool = multiprocessing.Pool(processes=num_processes)
    result_iter = pool.imap_unordered(run_classification_r_script, tasks)
    errors = False
    for retcode, output_prefix in result_iter:        
        if retcode == 0:
            result_table_file = output_prefix + ".classify.txt"
            shutil.copyfileobj(open(result_table_file), result_fh)
        else:
            errors = True
    pool.close()
    pool.join()
    result_fh.close()
    if errors:
        logging.error("Errors occurred during classification")

def sort_classification_results(input_file, output_file, tmp_dir):
    # sort classification results
    def sort_by_chrom_start(line):
        fields = line.strip().split('\t', 2)
        if fields[0] == "chrom":
            return chr(0), 0
        return fields[0], int(fields[1])
    batch_sort(input=input_file,
               output=output_file,
               key=sort_by_chrom_start,
               buffer_size=(1 << 21),
               tempdirs=[tmp_dir])

def classify_transcripts(classify_dir, num_processors, tmp_dir):
    # setup input and output files
    lib_counts_file = os.path.join(classify_dir, LIB_COUNTS_FILE)
    lib_counts_list = list(LibCounts.from_file(lib_counts_file))
    library_ids = [x.library_id for x in lib_counts_list]
    category_info_dict = {}
    for category_key in CATEGORIES:
        category_str = category_int_to_str[category_key]
        cinfo = CategoryInfo.create(library_ids, category_key, 
                                    category_str, classify_dir)
        category_info_dict[category_key] = cinfo
        # write input files for classifier
        logging.info("Writing classification input files category='%s'" % (cinfo.category_str))
        for transcripts in parse_gtf(open(cinfo.output_gtf_file)):
            for t in transcripts:
                library_id = t.attrs[GTFAttr.LIBRARY_ID]
                fields = get_classification_fields(t)
                # lookup file handle and open new file if necessary
                if not library_id in cinfo.result_fh_dict:
                    cinfo.result_fh_dict[library_id] = open(cinfo.result_file_dict[library_id], "w")        
                    print >>cinfo.result_fh_dict[library_id], '\t'.join(get_classification_header())
                # write to file
                print >>cinfo.result_fh_dict[library_id], '\t'.join(map(str, fields))        
        # close open file handles
        for fh in cinfo.result_fh_dict.itervalues():
            fh.close()
    for category_key, cinfo in category_info_dict.iteritems():
        classify_tasks = []
        for lib_counts in lib_counts_list:
            # see if can run classifier on this file
            if lib_counts.category_counts[category_key] > 0:
                filename = cinfo.result_file_dict[lib_counts.library_id]
                classify_tasks.append((lib_counts.library_id, filename))
        # run classification
        logging.info("Classifying transcripts category='%s'" % (cinfo.category_str))
        classify_category(cinfo, classify_tasks, num_processors, tmp_dir)
        # sort results
        sort_classification_results(cinfo.ctree_file, cinfo.sorted_ctree_file, tmp_dir)
        os.remove(cinfo.ctree_file)

def main():
    multiprocessing.freeze_support()
    # parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("-v", "--verbose", action="store_true", 
                        dest="verbose", default=False)
    parser.add_argument("-p", "--num-processors", type=int, 
                        dest="num_processors", default=1)
    parser.add_argument("classify_dir")
    args = parser.parse_args()
    # check command line parameters
    if not os.path.exists(args.classify_dir):
        parser.error("Classification directory %s not found" % (args.classify_dir))
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
    logging.info("verbose logging:         %s" % (args.verbose))
    logging.info("num processors:          %s" % (args.num_processors))
    logging.info("classify directory:      %s" % (args.classify_dir))
    tmp_dir = os.path.join(args.classify_dir, "tmp")
    if not os.path.exists(tmp_dir):
        logging.info("Creating tmp directory '%s'" % (tmp_dir))
        os.makedirs(tmp_dir)
    # run classification procedure        
    classify_transcripts(args.classify_dir, args.num_processors, tmp_dir)
    # cleanup
    if os.path.exists(tmp_dir):
        shutil.rmtree(tmp_dir)

if __name__ == "__main__":
    sys.exit(main())
