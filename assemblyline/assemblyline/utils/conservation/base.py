'''
Created on Jan 8, 2014

@author: mkiyer
'''
import os
import logging
import itertools
import glob
import re
import subprocess
import numpy as np

BIGWIG_TO_BEDGRAPH_BIN = 'bigWigToBedGraph'

def which(program):
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)
    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file
    return None

class BEDFeature(object):    
    __slots__ = ('chrom', 'tx_start', 'tx_end', 'name', 'score', 'strand',
                 'cds_start', 'cds_end', 'exon_count', 'block_starts', 
                 'block_sizes', 'exons', 'attr_fields')

    def __str__(self):
        fields = [self.chrom,
                  str(self.tx_start),
                  str(self.tx_end),
                  self.name,
                  str(self.score),
                  self.strand,
                  str(self.cds_start),
                  str(self.cds_end),
                  '0',
                  str(self.exon_count),
                  ','.join(map(str, self.block_sizes)) + ',',
                  ','.join(map(str, self.block_starts)) + ',']
        return '\t'.join(fields)

    @staticmethod
    def from_string(line):        
        if line is None:
            return None
        line = line.strip()
        if line.startswith('#'):
            logging.debug("skipping comment line: %s" % (line))
            return None
        if line.startswith('track'):
            logging.debug("skipping track header line: %s"  % (line))
            return None
        fields = line.split('\t')
        # first six fields are required
        g = BEDFeature()
        g.chrom = fields[0]
        g.tx_start = int(fields[1])
        g.tx_end = int(fields[2])
        if len(fields) <= 3:
            g.name = '%s:%d-%d' % (g.chrom, g.tx_start, g.tx_end)
        else:
            g.name = fields[3]
        if len(fields) <= 4:
            g.score = 0
        else:
            g.score = fields[4]
        if len(fields) <= 5:
            g.strand = '.'
        else:
            g.strand = fields[5]        
        if len(fields) <= 6:
            g.cds_start = g.tx_start
            g.cds_end = g.tx_end
            g.exon_count = 1
            g.exons = [(g.tx_start, g.tx_end)]
        else:
            g.cds_start = int(fields[6])
            g.cds_end = int(fields[7])
            g.exon_count = int(fields[9])
            g.block_sizes = map(int, fields[10].split(',')[:-1])
            g.block_starts = map(int, fields[11].split(',')[:-1])            
            g.exons = []
            for start, size in itertools.izip(g.block_starts, g.block_sizes):
                g.exons.append((g.tx_start + start, g.tx_start + start + size))
        if len(fields) <= 12:
            g.attr_fields = []
        else:
            g.attr_fields = fields[12:]
        return g

    @staticmethod
    def parse(line_iter):
        for line in line_iter:
            if not line:
                continue
            if not line.strip():
                continue            
            if line.startswith("#"):
                continue
            if line.startswith("track"):
                continue
            yield BEDFeature.from_string(line)
            
def find_bigwig_files(root_path, pattern):            
    # find bigwig files
    glob_pattern = pattern.replace("{{CHROM}}", "*")
    re_pattern = pattern.replace("{{CHROM}}", r'(.+)')
    chrom_bigwig_dict = {}
    for filename in glob.glob(os.path.join(root_path, glob_pattern)):
        m = re.match(re_pattern, filename)
        if m is not None:
            chrom = os.path.basename(m.group(1))
            chrom_bigwig_dict[chrom] = filename
    return chrom_bigwig_dict

def extract_bigwig_data(feature, bigwig_file):
    exon_arrays = []
    for exon_start,exon_end in feature.exons:
        exon_size = exon_end - exon_start
        arr = np.empty(exon_size, dtype=np.float)
        arr[:] = np.nan
        args = [BIGWIG_TO_BEDGRAPH_BIN, 
                '-chrom=%s' % (feature.chrom),
                '-start=%d' % (exon_start),
                '-end=%d' % (exon_end),
                bigwig_file,
                '/dev/stdout']
        p = subprocess.Popen(args, stdout=subprocess.PIPE)
        res = p.communicate()[0]
        for line in res.split('\n'):
            if not line:
                continue
            fields = line.strip().split('\t')
            line_start = int(fields[1])
            line_end = int(fields[2])
            value = float(fields[3])
            arr_start = line_start - exon_start
            arr_end = line_end - exon_start
            arr[arr_start:arr_end] = value
        if feature.strand == '-':
            arr = arr[::-1]
        exon_arrays.append(arr)
    if feature.strand == '-':
        exon_arrays.reverse()
    # compute summary statistics
    arr = np.concatenate(exon_arrays)
    return arr
