'''
Created on Dec 18, 2010

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
import itertools

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
        g.name = fields[3]
        if len(fields) <= 4:
            g.score = 0
            g.strand = '.'
        else:
            g.score = fields[4]
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
