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
import logging

# float precision threshold
FLOAT_PRECISION = 1e-10

# assemblyline GTF attributes
class GTFAttr:
    GENE_ID = 'gene_id'
    TRANSCRIPT_ID = 'transcript_id'
    SAMPLE_ID = 'sid'
    LIBRARY_ID = 'lid'
    REF = 'ref'
    TEST = 'tst'
    PCTRANK = 'pct'
    CATEGORY = 'cat'
    ANN_REF_ID = 'aid'
    ANN_COV_RATIO = 'acr'
    ANN_INTRON_RATIO = 'air'
    SCORE = 'score'
    MEAN_SCORE = 'avgscore'
    MEAN_PCTRANK = 'avgpct'
    MEAN_RECURRENCE = 'avgrecur'
    RESOLVED_STRAND = 'resolvedstrand'
    LOG10LR = 'log10lr'

class Category(object):
    # constant transcript category values
    SAME_STRAND = 0    
    OPP_STRAND = 1
    INTRONIC_SAME_STRAND = 2
    INTRONIC_OPP_STRAND = 3
    INTRONIC_AMBIGUOUS = 4
    INTERLEAVING = 5
    INTERGENIC = 6
    NUM_CATEGORIES = 7
    CATEGORIES = range(0, NUM_CATEGORIES)

    # groups of categories
    INTRONIC_LIKE = set([INTRONIC_SAME_STRAND, INTRONIC_AMBIGUOUS])
    INTERGENIC_LIKE = set([OPP_STRAND, INTRONIC_OPP_STRAND, INTERLEAVING, 
                           INTERGENIC])

    INT_TO_STR_DICT = {SAME_STRAND: 'same_strand',
                       OPP_STRAND: 'opp_strand',
                       INTRONIC_SAME_STRAND: 'intronic_same_strand',
                       INTRONIC_OPP_STRAND: 'intronic_opp_strand',
                       INTRONIC_AMBIGUOUS: 'intronic_ambiguous',
                       INTERLEAVING: 'interleaving',
                       INTERGENIC: 'intergenic'}
    
    @staticmethod
    def to_str(catint):
        return Category.INT_TO_STR_DICT[catint]

class CategoryCounts(object):
    def __init__(self):
        self.library_id = None
        self.counts = [0] * Category.NUM_CATEGORIES
    def to_fields(self):
        return [self.library_id] + self.counts
    @property
    def num_transcripts(self):
        return sum(self.counts)        
    @staticmethod
    def header_fields():
        fields = ["library_id"]
        fields.extend([Category.to_str(k) 
                       for k in range(0, Category.NUM_CATEGORIES)])
        return fields
    @staticmethod
    def from_line(line):
        fields = line.strip().split('\t')
        c = CategoryCounts()
        c.library_id = fields[0]
        c.counts = map(int, fields[1:])
        return c
    @staticmethod
    def from_file(filename):
        f = open(filename)
        f.next()
        for line in f:
            yield CategoryCounts.from_line(line)
        f.close()

class Library(object):
    fields = ('sample_id', 'library_id', 'gtf_file', 'bam_file')

    def __init__(self):
        for f in Library.fields:
            setattr(self, f, None)

    @staticmethod
    def from_fields(fields, field_dict=None):
        if field_dict is None:
            field_dict = dict((x,i) for i,x in enumerate(Library.fields))
        tbl = Library()
        for attrname in field_dict:
            setattr(tbl, attrname, fields[field_dict[attrname]])
        return tbl

    @staticmethod
    def from_file(filename):
        fh = open(filename)
        # header
        field_names = fh.next().strip().split('\t')
        field_dict = dict((x,i) for i,x in enumerate(field_names))
        # table rows
        for line in fh:
            fields = line.strip().split('\t')
            yield Library.from_fields(fields, field_dict)
        fh.close()
        
    def is_valid(self):
        if self.gtf_file is None:
            logging.error("Library %s GTF file is None" % (self.library_id))
            return False
        if not os.path.exists(self.gtf_file):
            logging.error("Library %s GTF file not found" % (self.library_id))
            return False
        if self.bam_file is None:
            logging.error("Library %s BAM file is None" % (self.library_id))
            return False
        if not os.path.exists(self.bam_file):
            logging.error("Library %s BAM file not found" % (self.library_id))
            return False
        return True

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

def cmp_strand(a,b):
    if a == "." or b == ".":
        return True
    return a == b