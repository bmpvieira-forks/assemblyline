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
import multiprocessing
import collections

import assemblyline
import assemblyline.lib.config as config
from assemblyline.lib.transcript import parse_gtf
from assemblyline.lib.base import CategoryCounts, Category, \
    GTFAttr, check_executable 
from assemblyline.lib.gtf import GTFFeature

# R script to call for classifying transcripts
_module_dir = assemblyline.__path__[0]
CLASSIFY_KDE2D_R_SCRIPT = os.path.join(_module_dir, "lib", 
                                       "classify_transcripts_kde2d.R")
DInfo = collections.namedtuple('DecisionInfo', ['pred', 'log10lr', 'is_test'])


def get_classify_header_fields():
    header = ["chrom", "start", "library_id", "t_id", "category", 
              "test", "ann_ref_id", "ann_cov_ratio", "ann_intron_ratio", 
              "mean_score", "mean_recurrence", "mean_pctrank", "pctrank",
              "length", "num_exons"]
    return header

def get_classify_fields(t):
    # setup list of annotation fields
    fields = [t.chrom,
              t.start,
              t.attrs[GTFAttr.LIBRARY_ID],
              t.attrs[GTFAttr.TRANSCRIPT_ID],
              t.attrs[GTFAttr.CATEGORY],
              t.attrs[GTFAttr.TEST],
              t.attrs[GTFAttr.ANN_REF_ID],
              t.attrs[GTFAttr.ANN_COV_RATIO],
              t.attrs[GTFAttr.ANN_INTRON_RATIO],
              t.attrs[GTFAttr.MEAN_SCORE],
              t.attrs[GTFAttr.MEAN_RECURRENCE],
              t.attrs[GTFAttr.MEAN_PCTRANK],
              t.attrs[GTFAttr.PCTRANK],
              t.length,
              len(t.exons)]
    return fields

def write_transcript_table(gtf_file, table_file):
    fileh = open(table_file, 'w')
    print >>fileh, '\t'.join(get_classify_header_fields())
    for transcripts in parse_gtf(open(gtf_file)):
        for t in transcripts:
            fields = get_classify_fields(t)
            print >>fileh, '\t'.join(map(str, fields))
    fileh.close()
    
def read_classify_stats(filename):
    field_dict = {}
    for line in open(filename):
        fields = line.strip().split('\t')
        field_dict[fields[0]] = fields[1:]
    return field_dict

def read_classify_decisions(filename, cutoff_type):
    fileh = open(filename)
    header_fields = fileh.next().strip().split('\t')
    t_id_ind = header_fields.index('t_id')
    category_ind = header_fields.index('category')
    test_ind = header_fields.index('test')
    num_exons_ind = header_fields.index('num_exons')
    log10lr_intronic_ind = header_fields.index('log10lr.intronic')
    log10lr_intergenic_ind = header_fields.index('log10lr.intergenic')
    if cutoff_type == "test":
        #intronic_pred_ind = header_fields.index("pred.test.intronic")
        intergenic_pred_ind = header_fields.index("pred.test.intergenic")
    else:
        #intronic_pred_ind = header_fields.index("pred.train.intronic")
        intergenic_pred_ind = header_fields.index("pred.train.intergenic")
    decision_dict = {}
    for line in fileh:
        fields = line.strip().split('\t')
        t_id = fields[t_id_ind]
        category = int(fields[category_ind])
        is_test = bool(int(fields[test_ind]))
        multi_exon = int(fields[num_exons_ind]) > 1        
        if is_test or category == Category.SAME_STRAND:
            intronic_log10lr = fields[log10lr_intronic_ind]
            intergenic_log10lr = fields[log10lr_intergenic_ind]
            if intronic_log10lr == "NA" and intergenic_log10lr == "NA":
                log10lr = "NA"
            elif intronic_log10lr != "NA" and intergenic_log10lr != "NA":
                log10lr = str(max(intronic_log10lr, intergenic_log10lr))
            elif intergenic_log10lr != "NA":
                log10lr = intergenic_log10lr
            else:
                log10lr = intronic_log10lr
            pred = True
        elif category in Category.INTRONIC_LIKE:
            pred = multi_exon
            log10lr = fields[log10lr_intronic_ind]
            #pred = fields[intronic_pred_ind] == "TRUE"
        elif category in Category.INTERGENIC_LIKE:
            pred = fields[intergenic_pred_ind] == "TRUE"
            log10lr = fields[log10lr_intergenic_ind]
        decision_dict[t_id] = DInfo(pred=pred,
                                    log10lr=log10lr,
                                    is_test=is_test)
    fileh.close()
    return decision_dict 

def classify_library_transcripts(args):
    library_id, output_dir, cutoff_type = args
    prefix = os.path.join(output_dir, library_id)
    # input files
    input_gtf_file = prefix + ".gtf"
    logfile = prefix + ".log"
    tablefile = prefix + '.inp.txt'
    # output files
    stats_file = prefix + ".stats.txt"
    output_res_file = prefix + ".out.txt"
    expr_gtf_file = prefix + ".expr.gtf"
    bkgd_gtf_file = prefix + ".bkgd.gtf"
    # write table of observations
    logging.debug("[STARTED]  library_id='%s'" % (library_id))
    write_transcript_table(input_gtf_file, tablefile)
    # run R script to do classification
    logfh = open(logfile, "w")
    retcode = subprocess.call(["Rscript", "--vanilla",
                               CLASSIFY_KDE2D_R_SCRIPT, 
                               prefix], 
                              stdout=logfh, stderr=logfh)
    logfh.close()
    if retcode != 0:
        logging.error("[FAILED]   library_id='%s'" % (library_id))
        return retcode, library_id
    # get library stats
    stats_field_dict = read_classify_stats(stats_file)
    has_tests = int(stats_field_dict["tests"][0]) > 0
    if not has_tests:
        cutoff_type = "train"
    # get transcript predictions
    decision_dict = read_classify_decisions(output_res_file, cutoff_type)
    # partition input into expressed vs background
    expr_fileh = open(expr_gtf_file, 'w')
    bkgd_fileh = open(bkgd_gtf_file, 'w')
    output_file_handles = [bkgd_fileh, expr_fileh]
    for feature in GTFFeature.parse(open(input_gtf_file)):
        t_id = feature.attrs[GTFAttr.TRANSCRIPT_ID]
        dinf = decision_dict[t_id]
        feature.attrs[GTFAttr.LOG10LR] = dinf.log10lr
        fileh = output_file_handles[int(dinf.pred)]
        print >>fileh, str(feature)
    for fileh in output_file_handles:
        fileh.close()
    logging.debug("[FINISHED] library_id='%s'" % (library_id))
    return retcode, library_id

def classify_transcripts(results, cutoff_type, num_processors):
    # read library category statistics
    counts_list = list(CategoryCounts.from_file(results.category_counts_file))
    # get tasks
    tasks = []
    for countsobj in counts_list:
        library_id = countsobj.library_id
        tasks.append((library_id, results.classify_dir, cutoff_type))
    # use multiprocessing to parallelize
    num_processes = max(1, num_processors - 1)
    pool = multiprocessing.Pool(processes=num_processes)
    result_iter = pool.imap_unordered(classify_library_transcripts, tasks)
    errors = False
    expressed_gtf_files = []
    background_gtf_files = []
    library_ids = []
    for retcode, library_id in result_iter:
        if retcode == 0:
            library_ids.append(library_id)
            prefix = os.path.join(results.classify_dir, library_id)
            expressed_gtf_files.append(prefix + ".expr.gtf")
            background_gtf_files.append(prefix + ".bkgd.gtf")
        else:
            errors = True
    pool.close()
    pool.join()
    if errors:
        logging.error("Errors occurred during classification")
    return int(errors)

def main():
    multiprocessing.freeze_support()
    # parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("-v", "--verbose", action="store_true", 
                        dest="verbose", default=False)
    parser.add_argument("-p", "--num-processors", type=int, 
                        dest="num_processors", default=1)
    parser.add_argument("--cutoff-type", dest="cutoff_type",
                        choices=["train", "test"], default="test",
                        help="Whether to choose classification "
                        "cutoff based on 'test' data or 'train' data "
                        "[default=%(default)s]")
    parser.add_argument("run_dir")
    args = parser.parse_args()
    # check command line parameters
    if not os.path.exists(args.run_dir):
        parser.error("Run directory %s not found" % (args.run_dir))
    # check command line parameters
    if not os.path.exists(CLASSIFY_KDE2D_R_SCRIPT):
        parser.error("Classification R script not found")
    if not check_executable("Rscript"):
        parser.error("Rscript binary not found")        
    num_processors = max(1, args.num_processors)
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
    logging.info("num processors:   %d" % (args.num_processors))
    logging.info("verbose logging:  %s" % (args.verbose))
    logging.info("cutoff type:      %s" % (args.cutoff_type))
    logging.info("----------------------------------")   
    # setup results
    results = config.AssemblylineResults(args.run_dir)
    return classify_transcripts(results, args.cutoff_type, 
                                num_processors)

if __name__ == "__main__":
    sys.exit(main())
