'''
Created on Feb 15, 2011

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
import os
import subprocess
import math
import collections

# constant file names
LIB_COUNTS_FILE = "lib_counts.txt"

# constant transcript category values
SENSE = 0
ANTISENSE = 1
INTRONIC = 2
INTERGENIC = 3
CATEGORIES = [SENSE, ANTISENSE, INTRONIC, INTERGENIC]
category_int_to_str = {SENSE: "sense",
                       ANTISENSE: "antisense",
                       INTRONIC: "intronic",
                       INTERGENIC: "intergenic"}
category_str_to_int = dict((v,k) for k,v in category_int_to_str.items())

# annotation GTF attributes
ANNOTATED = 'ann'
CATEGORY = 'cat'
MEAN_SCORE = 'mean_score'
MEAN_RECURRENCE = 'mean_recurrence'

# decision codes
ANN_EXPR = 0
ANN_BKGD = 1
UNANN_EXPR = 2
UNANN_BKGD = 3
SKIPPED = 4

class CategoryInfo():
    def __init__(self):
        self.category_key = None
        self.category_str = None
        self.library_ids = set()
        self.output_dir = None
        self.ctree_dir = None
        self.result_file_dict = {}
        self.result_fh_dict = {}
        self.cutoff_file_dict = {}
        self.output_gtf_file = None
        self.output_gtf_fh = None
        self.ctree_file = None
        self.sorted_ctree_file = None
        self.ann_expr_gtf_file = None
        self.unann_expr_gtf_file = None
        self.ann_bkgd_gtf_file = None
        self.unann_bkgd_gtf_file = None
        self.skipped_gtf_file = None
        self.decision_file_dict = {}
        self.decision_fh_dict = {}
        self.pred_stats_file = None

    @staticmethod
    def create(library_ids, category_key, category_str, output_dir):
        # make category info object
        cinfo = CategoryInfo()
        cinfo.category_key = category_key
        cinfo.category_str = category_str
        cinfo.library_ids = set(library_ids)
        category_dir = os.path.join(output_dir, category_str)
        if not os.path.exists(category_dir):
            os.makedirs(category_dir)
        cinfo.output_dir = category_dir
        ctree_dir = os.path.join(category_dir, "ctrees")
        if not os.path.exists(ctree_dir):
            os.makedirs(ctree_dir)
        cinfo.ctree_dir = ctree_dir
        result_file_dict = collections.defaultdict(lambda: {})
        cutoff_file_dict = collections.defaultdict(lambda: {})
        for library_id in library_ids:
            filename = os.path.join(ctree_dir, "%s.txt" % (library_id))
            result_file_dict[library_id] = filename
            filename = os.path.join(ctree_dir, "%s.cutoffs.txt" % (library_id))
            cutoff_file_dict[library_id] = filename
        cinfo.result_file_dict = result_file_dict
        cinfo.cutoff_file_dict = cutoff_file_dict
        cinfo.output_gtf_file = os.path.join(category_dir, "transcripts.gtf")
        cinfo.ctree_file = os.path.join(category_dir, "transcripts.classify.txt")
        cinfo.sorted_ctree_file = os.path.join(category_dir, "transcripts.classify.srt.txt")
        # filtered result files
        cinfo.ann_expr_gtf_file = os.path.join(category_dir, "ann_expr.gtf")
        cinfo.unann_expr_gtf_file = os.path.join(category_dir, "unann_expr.gtf")
        cinfo.ann_bkgd_gtf_file = os.path.join(category_dir, "ann_bkgd.gtf")
        cinfo.unann_bkgd_gtf_file = os.path.join(category_dir, "unann_bkgd.gtf")
        cinfo.skipped_gtf_file = os.path.join(category_dir, "skipped.gtf")
        cinfo.decision_file_dict = {ANN_EXPR: cinfo.ann_expr_gtf_file,
                                    ANN_BKGD: cinfo.ann_bkgd_gtf_file,
                                    UNANN_EXPR: cinfo.unann_expr_gtf_file,
                                    UNANN_BKGD: cinfo.unann_bkgd_gtf_file,
                                    SKIPPED: cinfo.skipped_gtf_file}
        cinfo.pred_stats_file = os.path.join(category_dir, "pred_stats.txt") 
        return cinfo

class LibCounts(object):
    def __init__(self):
        self.library_id = None
        self.annotated_counts = [0, 0]
        self.category_counts = [0] * len(CATEGORIES)
    @staticmethod
    def header_fields():
        fields = ["library_id", "unannotated", "annotated"]
        fields.extend([category_int_to_str[k] for k in CATEGORIES])
        return fields
    def to_fields(self):
        return [self.library_id] + self.annotated_counts + self.category_counts
    @staticmethod
    def from_line(line):
        fields = line.strip().split('\t')
        c = LibCounts()
        c.library_id = fields[0]
        c.annotated_counts = map(int, fields[1:3])
        c.category_counts = map(int, fields[3:])
        return c
    @staticmethod
    def from_file(filename):
        f = open(filename)
        f.next()
        for line in f:
            yield LibCounts.from_line(line)

def get_classification_header():
    header = ["chrom", "start", "library_id", "t_id", "annotated", 
              "category", "length", "num_exons", "score", 
              "mean_score", "mean_recurrence"]
    return header

def get_classification_result_header():
    header = get_classification_header()
    header.append("pred")
    return header

imax2 = lambda x,y: x if x>=y else y
imin2 = lambda x,y: x if x<=y else y

def float_check_nan(x):
    x = float(x)
    if math.isnan(x):
        return 0.0
    return x

def check_executable(filename):
    # check that samtools binary exists
    devnullfh = open(os.devnull, 'w')        
    try:
        subprocess.call([filename], stdout=devnullfh, stderr=devnullfh)
    except OSError:
        return False
    devnullfh.close()
    return True

