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
from assemblyline.lib.base import CategoryStats, Category, \
    GTFAttr, check_executable, FileHandleCache, BufferedFileSplitter
from assemblyline.lib.gtf import GTFFeature, merge_sort_gtf_files

# R script to call for classifying transcripts
_module_dir = assemblyline.__path__[0]
CLASSIFY_R_SCRIPT = os.path.join(_module_dir, "lib", "classify_transcripts.R")
DInfo = collections.namedtuple('DecisionInfo', ['pred', 'log10lr', 'is_test'])


def get_classify_header_fields():
    header = ["chrom", "start", "library_id", "t_id", "category", 
              "test", "ann_ref_id", "ann_cov_ratio", "ann_intron_ratio", 
              "mean_recurrence", "score", "pctrank", "length", "num_exons"]
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
              t.attrs[GTFAttr.MEAN_RECURRENCE],
              t.attrs[GTFAttr.SCORE],
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
    
def read_classify_info(filename):
    field_dict = {}
    for line in open(filename):
        fields = line.strip().split('\t')
        field_dict[fields[0]] = fields[1:]
    return field_dict

def read_classify_decisions(filename):
    fileh = open(filename)
    header_fields = fileh.next().strip().split('\t')
    t_id_ind = header_fields.index('t_id')
    category_ind = header_fields.index('category')
    test_ind = header_fields.index('test')
    num_exons_ind = header_fields.index('num_exons')
    log10lr_intronic_ind = header_fields.index('log10lr.intronic')
    log10lr_intergenic_ind = header_fields.index('log10lr.intergenic')
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
            if multi_exon:
                pred = True
            else:
                # For now do not allow single exon intronic same-stranded 
                # transcripts
                pred = False
                #pred = fields[intronic_pred_ind] == "TRUE"
            log10lr = fields[log10lr_intronic_ind]
        elif category in Category.INTERGENIC_LIKE:
            pred = multi_exon
            if multi_exon:
                pred = True
            else:
                pred = fields[intergenic_pred_ind] == "TRUE"
            log10lr = fields[log10lr_intergenic_ind]
        decision_dict[t_id] = DInfo(pred=pred,
                                    log10lr=log10lr,
                                    is_test=is_test)
    fileh.close()
    return decision_dict 

def classify_library_transcripts(args):
    library_id, output_dir = args
    prefix = os.path.join(output_dir, library_id)
    # input files
    input_gtf_file = prefix + ".gtf"
    logfile = prefix + ".log"
    tablefile = prefix + '.inp.txt'
    # output files
    info_file = prefix + ".info.txt"
    output_res_file = prefix + ".out.txt"
    expr_gtf_file = prefix + ".expr.gtf"
    bkgd_gtf_file = prefix + ".bkgd.gtf"
    # write table of observations
    logging.debug("[STARTED]  library_id='%s'" % (library_id))
    write_transcript_table(input_gtf_file, tablefile)
    # run R script to do classification
    logfh = open(logfile, "w")
    retcode = subprocess.call(["Rscript", "--vanilla",
                               CLASSIFY_R_SCRIPT, 
                               prefix], 
                              stdout=logfh, stderr=logfh)
    logfh.close()
    if retcode != 0:
        logging.error("[FAILED]   library_id='%s'" % (library_id))
        return retcode, library_id
    # get library stats
    #info_field_dict = read_classify_info(info_file)
    #has_tests = int(info_field_dict["tests"][0]) > 0
    # get transcript predictions
    decision_dict = read_classify_decisions(output_res_file)
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

def classify_transcripts(results, num_processors):
    # read library category statistics
    stats_list = list(CategoryStats.from_file(results.category_stats_file))
    # get tasks
    tasks = []
    for statsobj in stats_list:
        library_id = statsobj.library_id
        tasks.append((library_id, results.classify_dir))
    # use multiprocessing to parallelize
    pool = multiprocessing.Pool(processes=num_processors)
    result_iter = pool.imap_unordered(classify_library_transcripts, tasks)
    errors = False
    library_ids = []
    for retcode, library_id in result_iter:
        if retcode == 0:
            library_ids.append(library_id)
        else:
            errors = True
    pool.close()
    pool.join()
    if errors:
        logging.error("Errors occurred during classification")
    return int(errors)

def merge_transcripts(results):
    # read library category statistics
    stats_list = list(CategoryStats.from_file(results.category_stats_file))
    library_ids = []
    expressed_gtf_files = []
    background_gtf_files = []
    for statsobj in stats_list:
        library_id = statsobj.library_id
        library_ids.append(library_id)
        prefix = os.path.join(results.classify_dir, library_id)
        expressed_gtf_files.append(prefix + ".expr.gtf")
        background_gtf_files.append(prefix + ".bkgd.gtf")
    library_id_map = {}
    for line in open(results.library_id_map):
        fields = line.strip().split('\t')
        library_id_map[fields[0]] = fields[1]
    # make a classification report
    logging.info("Writing classification report")
    fileh = open(results.classify_report_file, 'w')
    header_fields = ["library_id", "library_name", "category", 
                     "train.auc", "test.auc", "train.cutoff", "train.tp",
                     "train.fp", "train.fn", "train.tn", "train.sens", "train.spec", 
                     "train.balacc", "test.tp", "test.fp", "test.fn", "test.tn", 
                     "test.sens", "test.spec", "test.balacc"]
    print >>fileh, '\t'.join(header_fields)
    for library_id in library_ids:
        prefix = os.path.join(results.classify_dir, library_id)
        library_name = library_id_map[library_id]
        intergenic_perf_file = prefix + ".intergenic.perf.txt"
        intronic_perf_file = prefix + ".intronic.perf.txt"
        input_fileh = open(intergenic_perf_file)
        input_fileh.next()
        for line in input_fileh:
            fields = ([library_id, library_name, "intergenic"] 
                      + line.strip().split('\t'))
            print >>fileh, '\t'.join(fields)
        input_fileh.close()
        input_fileh = open(intronic_perf_file)
        input_fileh.next()        
        for line in input_fileh:
            fields = ([library_id, library_name, "intronic"] 
                      + line.strip().split('\t'))
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
    return 0

#def split_gtf_file(gtf_file, split_dir, ref_gtf_file, category_stats_file,
#                   bufsize=(1 << 30)):
#    # split input gtf by library and mark test ids
#    keyfunc = lambda myid: os.path.join(split_dir, "%s.gtf" % (myid))
#    cache = FileHandleCache(keyfunc)
#    ref_fileh = open(ref_gtf_file, 'w')
#    stats_dict = collections.defaultdict(lambda: CategoryStats())
#    logging.info("Splitting transcripts by library")
#    for f in GTFFeature.parse(open(gtf_file)):
#        is_ref = bool(int(f.attrs[GTFAttr.REF]))
#        if is_ref:
#            print >>ref_fileh, str(f)
#            continue
#        library_id = f.attrs[GTFAttr.LIBRARY_ID]
#        # keep statistics
#        if f.feature_type == 'transcript':
#            category = int(f.attrs[GTFAttr.CATEGORY])
#            score = float(f.attrs[GTFAttr.SCORE])         
#            statsobj = stats_dict[library_id]
#            statsobj.library_id = library_id
#            statsobj.counts[category] += 1
#            statsobj.signal[category] += score
#        # write features from each library to separate files
#        fileh = cache.get_file_handle(library_id)
#        print >>fileh, str(f)
#    # close open file handles
#    ref_fileh.close()
#    cache.close()
#    logging.debug("File handle cache hits: %d" % (cache.hits))
#    logging.debug("File handle cache misses: %d" % (cache.misses))
#    # write library category statistics
#    logging.info("Writing category statistics")
#    fh = open(category_stats_file, "w")
#    print >>fh, '\t'.join(CategoryStats.header_fields())
#    for statsobj in stats_dict.itervalues():
#        fields = statsobj.to_fields()
#        print >>fh, '\t'.join(map(str, fields))
#    fh.close()

def split_gtf_file(gtf_file, split_dir, ref_gtf_file, category_stats_file,
                   bufsize=(1 << 30)):
    # split input gtf by library and mark test ids
    keyfunc = lambda myid: os.path.join(split_dir, "%s.gtf" % (myid))
    bufobj = BufferedFileSplitter(keyfunc, bufsize)
    ref_fileh = open(ref_gtf_file, 'w')
    stats_dict = collections.defaultdict(lambda: CategoryStats())
    logging.info("Splitting transcripts by library")
    for line in open(gtf_file):
        f = GTFFeature.from_string(line)
        is_ref = bool(int(f.attrs[GTFAttr.REF]))
        if is_ref:
            print >>ref_fileh, str(f)
            continue
        library_id = f.attrs[GTFAttr.LIBRARY_ID]
        # keep statistics
        if f.feature_type == 'transcript':
            is_test = bool(int(f.attrs[GTFAttr.TEST]))
            if is_test:
                category = Category.SAME_STRAND
            else:
                category = int(f.attrs[GTFAttr.CATEGORY])
            score = float(f.attrs[GTFAttr.SCORE])         
            statsobj = stats_dict[library_id]
            statsobj.library_id = library_id
            statsobj.counts[category] += 1
            statsobj.signal[category] += score
        # write features from each library to separate files
        bufobj.write(library_id, line)
    # close open file handles
    ref_fileh.close()
    bufobj.close()
    logging.debug("Buffer flushes: %d" % (bufobj.flushes))
    # write library category statistics
    logging.info("Writing category statistics")
    fh = open(category_stats_file, "w")
    print >>fh, '\t'.join(CategoryStats.header_fields())
    for statsobj in stats_dict.itervalues():
        fields = statsobj.to_fields()
        print >>fh, '\t'.join(map(str, fields))
    fh.close()

def main():
    multiprocessing.freeze_support()
    # parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("-v", "--verbose", action="store_true", 
                        dest="verbose", default=False)
    parser.add_argument("--bufsize", dest="bufsize", type=int,
                        default=(1 << 30),
                        help="Size of buffer when splitting GTF file")
    parser.add_argument("-p", "--num-processors", type=int, 
                        dest="num_processors", default=1)
    parser.add_argument("run_dir")
    args = parser.parse_args()
    # check command line parameters
    if not os.path.exists(args.run_dir):
        parser.error("Run directory %s not found" % (args.run_dir))
    # check command line parameters
    if not os.path.exists(CLASSIFY_R_SCRIPT):
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
    logging.info("buffer size:      %d" % (args.bufsize))
    logging.info("verbose logging:  %s" % (args.verbose))
    logging.info("----------------------------------")   
    # setup results
    results = config.AssemblylineResults(args.run_dir)
    if not os.path.exists(results.classify_dir):
        os.makedirs(results.classify_dir)
    # split gtf file
    split_gtf_file(results.annotated_transcripts_gtf_file, 
                   results.classify_dir,
                   results.ref_gtf_file,
                   results.category_stats_file,
                   args.bufsize)
    # run classification
    retcode = classify_transcripts(results, num_processors)
    if retcode != 0:
        logging.error("ERROR")
        return retcode
    # merge results
    retcode = merge_transcripts(results)
    logging.info("Done")
    return retcode

if __name__ == "__main__":
    sys.exit(main())
