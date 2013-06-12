'''
Created on Jan 19, 2013

@author: mkiyer
'''

class BEDFeature(object):
    __slots__ = ('chrom', 'tx_start', 'tx_end', 'name', 'score', 'strand',
                 'cds_start', 'cds_end', 'exon_count', 'exons', 'introns')
    
    def __init__(self):
        self.chrom = None
        self.tx_start = 0
        self.tx_end = 0
        self.name = None
        self.score = 0.0
        self.strand = "."
        self.cds_start = 0
        self.cds_end = 0
        self.exon_count = 0
        self.exons = []
        self.introns = []

    @staticmethod
    def parse_line(line):
        if line is None:
            return None
        line = line.strip()
        if line.startswith('#'):
            return None
        if line.startswith('track'):
            return None
        fields = line.split('\t')
        f = BEDFeature()
        f.chrom = fields[0]
        f.tx_start = int(fields[1])
        f.tx_end = int(fields[2])
        f.name = fields[3]
        if len(fields) < 5:
            f.exon_count = 1
            f.exons = [(f.tx_start,f.tx_end)]
            return f
        f.score = fields[4]
        f.strand = fields[5]
        if len(fields) < 7:
            f.exon_count = 1
            f.exons = [(f.tx_start,f.tx_end)]
            return f
        f.cds_start = int(fields[6])
        f.cds_end = int(fields[7])
        if len(fields) < 10:
            f.exon_count = 1
            f.exons = [(f.tx_start,f.tx_end)]
            return f
        f.exon_count = int(fields[9])
        block_sizes = map(int, fields[10].split(',')[:-1])
        block_starts = map(int, fields[11].split(',')[:-1])        
        exon_starts = [(f.tx_start + start) for start in block_starts]        
        exon_ends = [(start + size) for start, size in zip(exon_starts, block_sizes)]
        f.exons = zip(exon_starts, exon_ends)
        f.introns = zip(exon_ends, exon_starts[1:])        
        return f
    
    @staticmethod
    def parse(line_iter):
        for line in line_iter:
            g = BEDFeature.parse_line(line)
            if g is None:
                continue
            yield g