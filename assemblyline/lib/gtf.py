'''
Created on Nov 2, 2010

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
import shutil

GTF_EMPTY_FIELD = '.'
GTF_ATTR_SEP = ';'
GTF_ATTR_TAGVALUE_SEP = ' '

class GTFAttr:
    GENE_ID = 'gene_id'
    TRANSCRIPT_ID = 'transcript_id'
    COHORT_ID = 'cohort_id'
    SAMPLE_ID = 'sample_id'
    LIBRARY_ID = 'library_id'
    REF = 'ref'

class GTFError(Exception):
    pass

def sort_gtf(filename, output_file, tmp_dir=None):
    args = ["sort"]
    if tmp_dir is not None:
        args.extend(["-T", tmp_dir])
    args.extend(["-k1,1", "-k4,4n", "-k3,3r", filename])
    myenv = os.environ.copy()
    myenv["LC_ALL"] = "C"
    return subprocess.call(args, stdout=open(output_file, "w"), env=myenv)

def merge_sort_gtf_files(gtf_files, output_file, tmp_dir=None):
    tmp_file = os.path.splitext(output_file)[0] + ".unsorted.gtf"
    outfh = open(tmp_file, "w")
    for filename in gtf_files:
        shutil.copyfileobj(open(filename), outfh)
    outfh.close()
    sort_gtf(tmp_file, output_file, tmp_dir)
    os.remove(tmp_file)

def parse_loci(line_iter):
    '''
    requires that GTF file has been sorted and formatted such that a
    single 'transcript' feature appears before individual 'exon' 
    features such that transcript boundaries can be ascertained. this
    greatly simplifies parsing. using this function without appropriately 
    formatted GTF files will result in undefined behavior
    '''
    def window_overlap(a, b):
        if a[0] != b[0]:
            return False
        return (a[1] <= b[2]) and (b[1] <= a[2])
    def get_intervals(line_iter):
        for line in line_iter:
            # read the essential part of the GTF line
            line = line.rstrip()
            fields = line.split('\t', 5)
            seqid = fields[0]
            start = int(fields[3])-1
            end = int(fields[4])
            yield seqid, start, end, line
    try:
        interval_iter = get_intervals(line_iter)
        # initialize window
        seqid, start, end, line = interval_iter.next()
        window = [line]
        window_range = (seqid, start, end)
        # separate into loci
        for seqid, start, end, line in interval_iter:
            # check if next transcript is outside current window
            interval = (seqid, start, end)
            if not window_overlap(interval, window_range):
                # yield current window
                yield window
                # reset window
                window = [line]
                window_range = (seqid, start, end)
            else:
                # add transcript to window
                window.append(line)
                newstart = (start if start < window_range[1] else window_range[1])
                newend = (end if end > window_range[2] else window_range[2])
                window_range = (seqid, newstart, newend)
    except StopIteration:
        pass
    # yield last window
    if len(window) > 0:
        yield window

class GTFFeature(object):
    '''
    1. seqname - The name of the sequence. Must be a chromosome or scaffold.
    2. source - The program that generated this feature.
    3. feature - The name of this type of feature. Some examples of standard feature types are "CDS", "start_codon", "stop_codon", and "exon".
    4. start - The starting position of the feature in the sequence. The first base is numbered 1.
    5. end - The ending position of the feature (inclusive).
    6. score - A score between 0 and 1000. If the track line useScore attribute is set to 1 for this annotation data set, the score value will determine the level of gray in which this feature is displayed (higher numbers = darker gray). If there is no score value, enter ".".
    7. strand - Valid entries include '+', '-', or '.' (for don't know/don't care).
    8. phase - If the feature is a coding exon, frame should be a number between 0-2 that represents the reading frame of the first base. If the feature is not a coding exon, the value should be '.'.

    chr1    Cufflinks       transcript      136546  137059  1000    .       .       gene_id "VCAP_SHEZH2.657699"; transcript_id "VCAP_SHEZH2.657699.1"; FPKM "100.7219943204"; frac "1.000000"; conf_lo "80.649925"; conf_hi "120.794064"; cov "2.198209";
    '''
    __slots__ = ('seqid', 'source', 'feature_type', 'start', 'end', 'score', 'strand', 'phase', 'attrs') 

    def __str__(self):
        line = [self.seqid,
                self.source,
                self.feature_type,
                # convert to 1-based intervals
                str(self.start + 1),
                str(self.end),
                str(self.score),
                str(self.strand),
                self.phase]
        attr_str = ' '.join('%s "%s";' % (k, v) for (k, v) in self.attrs.iteritems())
        line.append(attr_str)
        return '\t'.join(line)
        
    @staticmethod
    def from_string(line, attr_defs=None):
        f = GTFFeature()        
        # read the GTF line
        fields = line.strip().split('\t')
        f.seqid = fields[0]
        f.source = fields[1]
        f.feature_type = fields[2]
        # convert from 1-based (inclusive) to 0-based (exclusive) intervals
        f.start = int(fields[3])-1
        f.end = int(fields[4])
        f.score = 0 if (fields[5] == '.') else float(fields[5])
        strand = fields[6]
        if not (strand == '+' or strand == '-'):
            strand = GTF_EMPTY_FIELD
        f.strand = strand        
        f.phase = fields[7]        
        attrs = {}
        if fields[8] != GTF_EMPTY_FIELD:
            attr_strings = fields[8].split(GTF_ATTR_SEP)
            for a in attr_strings:
                a = a.strip()
                if len(a) == 0:
                    continue
                tag, value = a.split(GTF_ATTR_TAGVALUE_SEP, 1)
                # remove quotes
                value = value.split('"')[1]
                # apply parsing function
                if (attr_defs != None) and (tag in attr_defs) and (attr_defs[tag] != None):
                    value = attr_defs[tag](value)
                attrs[tag] = value
        f.attrs = attrs
        return f

    @staticmethod
    def parse(line_iter, attr_defs=None):
        for line in line_iter:
            # read the GTF line
            if not line:
                continue
            if not line.strip():
                continue            
            if line.startswith("#"):
                continue
            yield GTFFeature.from_string(line, attr_defs)
