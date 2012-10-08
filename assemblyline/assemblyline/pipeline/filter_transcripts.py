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
import heapq
import collections
import shutil

import assemblyline
from assemblyline.lib.gtf import GTFFeature, sort_gtf
from assemblyline.lib.base import GTFAttr

from classify_transcripts import CategoryInfo, LibCounts, \
    LIB_COUNTS_FILE, category_int_to_str, CATEGORIES, \
    get_classification_result_header

# default parameters
DEFAULT_MIN_PREC = 0.0
DEFAULT_MIN_REC = 0.0
DEFAULT_MIN_SPEC = 0.5

# decision codes
ANN_EXPR = 0
ANN_BKGD = 1
UNANN_EXPR = 2
UNANN_BKGD = 3
SKIPPED = 4

# output files
EXPR_GTF_FILE = "expressed.gtf"
BKGD_GTF_FILE = "background.gtf"
ANN_GTF_FILE = "annotated.gtf"

# GTF attributes to add to final output files
GTF_ATTRS_TO_RETAIN = ('annotated', 'category', 
                       'mean_score', 'mean_recurrence', 'pred')

class ClassificationResult(object):
    __slots__ = get_classification_result_header()
    @staticmethod
    def from_line(line):
        fields = line.strip().split('\t')
        if len(fields) != 12:
            print fields
        c = ClassificationResult()
        c.chrom = fields[0]
        c.start = int(fields[1])
        c.library_id = fields[2]
        c.t_id = fields[3]
        c.annotated = int(fields[4])
        c.category = int(fields[5])
        c.length = int(fields[6])
        c.num_exons = int(fields[7])
        c.score = float(fields[8])
        c.mean_score = float(fields[9])
        c.mean_recurrence = float(fields[10])
        c.pred = float(fields[11])
        return c

def choose_cutoff(cutoff_table_file, min_prec, min_rec, min_spec, 
                  opt_variable):
    f = open(cutoff_table_file)
    # skip header
    header_fields = f.next().strip().split('\t')
    cutoff_ind = header_fields.index("cutoff")
    prec_ind = header_fields.index("prec")
    rec_ind = header_fields.index("rec")
    spec_ind = header_fields.index("spec")
    F_ind = header_fields.index("F")
    balacc_ind = header_fields.index("balacc")
    opt_ind = header_fields.index(opt_variable)
    # read entire table
    lines = [map(float, line.strip().split('\t')) for line in f]
    f.close()
    # find cutoffs that pass minimum thresholds
    valid_lines = []
    for i,fields in enumerate(lines):
        prec = fields[prec_ind]
        rec = fields[rec_ind]
        spec = fields[spec_ind]
        if ((prec >= min_prec) and (rec >= min_rec) and 
            (spec >= min_spec)):
            valid_lines.append(i)
    if len(valid_lines) == 0:
        return None, None, None, None, None, None
    # find cutoff with maximum F measure among valid cutoffs
    best_ind = None
    best_val = 0.0
    for i in valid_lines:
        fields = lines[i]
        val = fields[opt_ind]
        if (best_ind is None) or (val >= best_val):
            best_ind = i
            best_val = val 
    fields = lines[best_ind]
    return (fields[cutoff_ind], fields[prec_ind], fields[rec_ind],
            fields[spec_ind], fields[F_ind], fields[balacc_ind])

def run_filter(cinfo, cutoff_dict):
    # maintain predictions in a dictionary where key is transcript id and 
    # value is boolean prediction decision for transcript.
    t_id_decisions = {}
    # maintain result objects in dictionary keyed by transcript id
    t_id_results = {}
    # maintain heap queues that keeps track of the last transcript position 
    # on each chromosome. prediction decisions only need to be remembered
    # until the parsing goes past the end of the transcript (all exons 
    # accounted for)
    decision_heapqs = collections.defaultdict(lambda: [])
    # maintain heapq that keeps track of transcript position of
    # result objects.  results only need to be remembered until parsing
    # goes past the chrom/start position of the transcript (all transcripts
    # accounted for)
    result_heapqs = collections.defaultdict(lambda: [])
    # read result file and gtf file in sync
    result_fh = open(cinfo.sorted_ctree_file)
    result_fh.next()
    gtf_fh = open(cinfo.output_gtf_file)
    # open output files
    for decision,filename in cinfo.decision_file_dict.iteritems():
        cinfo.decision_fh_dict[decision] = open(filename, "w")
    # keep track of prediction statistics
    decision_stats = collections.defaultdict(lambda: 0)
    for feature in GTFFeature.parse(gtf_fh):
        # get transcript id used to lookup expressed/background 
        # prediction decision
        t_id = feature.attrs[GTFAttr.TRANSCRIPT_ID]
        # check top of the decision heapq and pop transcripts when parsing 
        # has gone past the end
        decision_heapq = decision_heapqs[feature.seqid]               
        while (len(decision_heapq) > 0) and (feature.start > decision_heapq[0][0]):
            smallest_end, smallest_t_id = heapq.heappop(decision_heapq)
            del t_id_decisions[smallest_t_id]
        # check top of result heapq and pop transcripts when parsing has gone
        # past the end
        result_heapq = result_heapqs[feature.seqid]
        while (len(result_heapq) > 0) and (feature.start > result_heapq[0][0]):
            result_start, result_t_id = heapq.heappop(result_heapq)
            del t_id_results[result_t_id]
        # parse transcript/exon features differently
        if feature.feature_type == "transcript":
            # parse results until this t_id is found (all results 
            # must stay valid until past this chrom/start location)
            while t_id not in t_id_results:
                result = ClassificationResult.from_line(result_fh.next())
                # add to heapq to remove results that are no longer useful
                heapq.heappush(result_heapqs[result.chrom], 
                               (result.start, result.t_id))
                # add to result dictionary
                t_id_results[result.t_id] = result
                # if current result position is beyond current transcript
                # position then we know that we are missing results for this
                # transcript and need to skip it
                if ((result.chrom != feature.seqid) or 
                    (result.start > feature.start)):
                    break
            if t_id not in t_id_results:
                #logging.warning("Skipping: library_id=%s t_id=%s "
                #                "chrom=%s start=%d " %
                #                (feature.attrs[GTFAttr.LIBRARY_ID], t_id,
                #                 feature.seqid, feature.start))
                decision = SKIPPED
            else:
                # lookup classification result and ensure that transcript_id 
                # attribute matches result id
                result = t_id_results[t_id]
                # lookup cutoff value for classification
                library_id = feature.attrs[GTFAttr.LIBRARY_ID]
                cutoff = cutoff_dict[library_id]
                feature.attrs["cutoff"] = cutoff
                is_expr = (result.pred >= cutoff)
                # retain certain results as transcript attributes
                for attr_name in GTF_ATTRS_TO_RETAIN:
                    feature.attrs[attr_name] = getattr(result, attr_name)
                # keep track of prediction decision and statistics
                # remember decision in dict so that it can be 
                # applied to the transcript exons as well
                if result.annotated:
                    if is_expr:
                        decision = ANN_EXPR
                    else:
                        decision = ANN_BKGD
                else:
                    if is_expr:
                        decision = UNANN_EXPR
                    else:
                        decision = UNANN_BKGD
            # push transcript end onto decision heap queue 
            # (decision must stay valid until past the end)
            heapq.heappush(decision_heapq, (feature.end, t_id))                                
            # keep track of decision to apply it to exon features
            t_id_decisions[t_id] = decision
        else:
            decision = t_id_decisions[t_id]
        # keep track of stats
        decision_stats[decision] += 1
        # output to separate files
        out_fh = cinfo.decision_fh_dict[decision]
        print >>out_fh, str(feature)
    # cleanup
    gtf_fh.close()
    result_fh.close()    
    for fh in cinfo.decision_fh_dict.itervalues():
        fh.close()
    return decision_stats

def filter_category(lib_counts_list, cinfo, min_prec, min_rec, min_spec, 
                    opt_variable, tmp_dir):
    logging.info("Choosing cutoffs category='%s'" % (cinfo.category_str))
    cutoff_dict = {}
    for lib_counts in lib_counts_list:
        # choose cutoff for each library
        library_id = lib_counts.library_id
        if lib_counts.category_counts[cinfo.category_key] == 0:
            logging.debug("\tlibrary_id='%s' category='%s' SKIPPED" % 
                          (library_id, cinfo.category_str))
            cutoff_dict[library_id] = 1.0
            continue
        cutoff_table_file = cinfo.cutoff_file_dict[library_id]
        cutoff, prec, rec, spec, F, balacc = \
            choose_cutoff(cutoff_table_file, min_prec, min_rec, 
                          min_spec, opt_variable)
        if cutoff is None:
            logging.debug("\tlibrary_id='%s' category='%s' SKIPPED" % 
                          (library_id, cinfo.category_str))
            cutoff_dict[library_id] = 1.0
        else:
            logging.debug("\tlibrary_id='%s' category='%s' cutoff=%f prec=%f "
                          "rec=%f spec=%f F=%f balacc=%f" % 
                          (library_id, cinfo.category_str, cutoff, prec, rec, 
                           spec, F, balacc))
            cutoff_dict[library_id] = cutoff
    # filter transcripts
    logging.info("Filtering transcripts category='%s'" % 
                 (cinfo.category_str))
    decision_stats = run_filter(cinfo, cutoff_dict)
    # compute statistics
    tp = decision_stats[ANN_EXPR]
    fp = decision_stats[UNANN_EXPR]
    fn = decision_stats[ANN_BKGD]
    tn = decision_stats[UNANN_BKGD]
    prec = tp / float(tp + fp)
    rec = tp / float(tp + fn)
    sens = rec
    spec = (tn / float(tn + fp))
    acc = (tp + tn) / float(tp + fp + tn + fn)
    balacc = (sens + spec) / 2.0
    if tp == 0:
        F = 0
    else:
        F = 2 * (prec * rec) / (prec + rec)    
    fh = open(cinfo.pred_stats_file, "w")
    print >>fh, "TP: %d" % (tp)
    print >>fh, "FN: %d" % (fn)
    print >>fh, "FP: %d" % (fp)
    print >>fh, "TN: %d" % (tn)
    print >>fh, "Skipped: %d" % (decision_stats[SKIPPED])
    print >>fh, "Precision: %f" % (prec)
    print >>fh, "Recall: %f" % (rec)
    print >>fh, "Sensitivity: %f" % (sens)
    print >>fh, "Specificity: %f" % (spec)
    print >>fh, "F: %f" % (F)
    print >>fh, "Accuracy: %f" % (acc)
    print >>fh, "Balanced accuracy: %f" % (balacc)
    fh.close()
    logging.info("\tTP=%d FN=%d FP=%d TN=%d" % (tp, fn, fp, tn))
    logging.info("\tSkipped: %d" % (decision_stats[SKIPPED]))
    logging.info("\tPrecision=%f Recall/Sensitivity=%f" % (prec, rec))
    logging.info("\tSpecificity=%f" % (spec))
    logging.info("\tF=%f" % F)
    logging.info("\tAccuracy=%f" % (acc))
    logging.info("\tBalanced accuracy=%f" % (balacc))
    return 0

def merge_sort_gtf_files(gtf_files, output_file):
    tmp_file = os.path.splitext(output_file)[0] + ".unsorted.gtf"
    outfh = open(tmp_file, "w")
    for filename in gtf_files:
        shutil.copyfileobj(open(filename), outfh)
    outfh.close()
    sort_gtf(tmp_file, output_file)
    os.remove(tmp_file)

def filter_transcripts(classify_dir, min_prec, min_rec, min_spec, 
                       opt_variable, tmp_dir):
    # setup input and output files
    lib_counts_file = os.path.join(classify_dir, LIB_COUNTS_FILE)
    lib_counts_list = list(LibCounts.from_file(lib_counts_file))
    # filter each category
    expr_gtf_files = []
    bkgd_gtf_files = []
    for category_key in CATEGORIES:
        category_str = category_int_to_str[category_key]
        category_lib_counts = [x for x in lib_counts_list
                               if x.category_counts[category_key] > 0]
        library_ids = [x.library_id for x in category_lib_counts]
        cinfo = CategoryInfo.create(library_ids, category_key, 
                                    category_str, classify_dir)
        filter_category(category_lib_counts, cinfo, min_prec, min_rec, 
                        min_spec, opt_variable, tmp_dir)
        expr_gtf_files.append(cinfo.unann_expr_gtf_file)
        bkgd_gtf_files.append(cinfo.unann_bkgd_gtf_file)
    # only need one set of annotated gtf files
    expr_gtf_files.extend([cinfo.ann_expr_gtf_file,
                           cinfo.ann_bkgd_gtf_file])
    # merge transcripts
    logging.info("Merging filtered transcripts")
    expr_gtf_file = os.path.join(classify_dir, EXPR_GTF_FILE)
    bkgd_gtf_file = os.path.join(classify_dir, BKGD_GTF_FILE)
    merge_sort_gtf_files(expr_gtf_files, expr_gtf_file)
    merge_sort_gtf_files(bkgd_gtf_files, bkgd_gtf_file)

def main():
    # parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("-v", "--verbose", action="store_true", 
                        dest="verbose", default=False)
    parser.add_argument("--prec", dest="min_prec", type=float,
                        default=DEFAULT_MIN_PREC)
    parser.add_argument("--rec", dest="min_rec", type=float,
                        default=DEFAULT_MIN_REC)
    parser.add_argument("--spec", dest="min_spec", type=float,
                        default=DEFAULT_MIN_SPEC)
    parser.add_argument("--optimize", dest="opt_variable",
                        choices=("F", "balacc"),
                        default="balacc")
    parser.add_argument("classify_dir")
    args = parser.parse_args()
    # check command line parameters
    if not os.path.exists(args.classify_dir):
        parser.error("Classify directory %s not found" % (args.classify_dir))
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
    logging.info("classify directory:    %s" % (args.classify_dir))
    logging.info("min precision:         %f" % (args.min_prec))
    logging.info("min recall:            %f" % (args.min_rec))
    logging.info("min specificity:       %f" % (args.min_spec))
    logging.info("optimize variable:     %s" % (args.opt_variable))
    # make tmp dir
    tmp_dir = os.path.join(args.classify_dir, "tmp")
    if not os.path.exists(tmp_dir):
        logging.info("Creating tmp directory '%s'" % (tmp_dir))
        os.makedirs(tmp_dir)
    filter_transcripts(args.classify_dir, args.min_prec, args.min_rec, 
                       args.min_spec, args.opt_variable, tmp_dir) 
    if os.path.exists(tmp_dir):
        shutil.rmtree(tmp_dir)
    return 0

if __name__ == '__main__':
    sys.exit(main())