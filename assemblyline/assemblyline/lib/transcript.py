'''
Created on Nov 30, 2011

@author: mkiyer
'''
from gtf import GTFFeature

# strand constants
POS_STRAND = 0
NEG_STRAND = 1
NO_STRAND = 2

# strand format conversion
_strand_int_to_str = ["+", "-", "."]
_strand_str_to_int = {"+": POS_STRAND,
                      "-": NEG_STRAND,
                      ".": NO_STRAND}
def strand_str_to_int(strand):
    return _strand_str_to_int[strand]
def strand_int_to_str(strand):
    return _strand_int_to_str[strand]

def interval_overlap(a, b):
    """return True if two intervals overlap, False otherwise"""
    return (a.start < b.end) and (b.start < a.end)

class Exon(object):
    __slots__ = ('start', 'end')
    def __init__(self, start, end):
        self.start = start
        self.end = end
    def __str__(self):
        return ("%s-%s" % (self.start, self.end))
    def __repr__(self):
        return ("<%s(start=%d end=%d)>" %
                (self.__class__.__name__, self.start, self.end))
    def __eq__(self, other):
        return (self.start == other.start) and (self.end == other.end)
    def __ne__(self, other):
        return (self.start != other.start) or (self.end != other.end)
    def __lt__(self, other):
        return self.start < other.start
    def __le__(self, other):        
        return self.start <= other.start
    def __gr__(self, other):
        return self.start > other.start
    def __ge__(self, other):
        return self.start >= other.start        
    def __hash__(self):
        return (self.start << 16) | (self.end)
    def is_overlapping(self, other):
        return interval_overlap(self, other)
    
class Transcript(object):
    __slots__ = ('chrom', 'start', 'end', 'strand', 'id', 'library', 
                 'sample', 'cohort', 'fpkm', 'cov', 'frac', 'exons', 
                 'attrs')

    def __init__(self):
        self.chrom = None
        self.start = -1
        self.end = -1
        self.strand = NO_STRAND
        self.id = None
        self.library = None
        self.sample = None
        self.cohort = None
        self.fpkm = 0.0
        self.cov = 0.0
        self.frac = 0.0
        self.exons = None
        self.attrs = None

    def __str__(self):
        return ("<%s(chrom='%s', start='%d', end='%d', strand='%s', "
                "id='%s', library='%s', cohort='%s', sample='%s', fpkm='%s', " 
                "cov='%s', frac='%s', length='%s', exons='%s'" % 
                (self.__class__.__name__, self.chrom, self.start, self.end, 
                 strand_int_to_str(self.strand), self.id, self.cohort, 
                 self.sample, self.fpkm, self.cov, self.frac, self.length, self.exons))

    @property
    def tx_start(self):
        return self.exons[0].start
    @property
    def tx_end(self):
        return self.exons[-1].end
    @property
    def length(self):
        return sum((e.end - e.start) for e in self.exons)
    @property
    def mass(self):
        return self.fpkm * float(self.length) * 1e-3
    @property
    def introns(self):
        return zip([e[1] for e in self.exons[:-1]],
                   [e[0] for e in self.exons[1:]])

    def to_bed12(self):
        block_sizes = []
        block_starts = []
        for e0, e1 in self.exons:
            block_starts.append(e0 - self.tx_start)
            block_sizes.append(e1 - e0)        
        # write
        s = '\t'.join([self.chrom, 
                       str(self.tx_start), 
                       str(self.tx_end),
                       str(self.id),
                       str(self.fpkm),
                       strand_int_to_str(self.strand),
                       str(self.tx_start),
                       str(self.tx_start),
                       '0',
                       str(len(self.exons)),
                       ','.join(map(str,block_sizes)) + ',',
                       ','.join(map(str,block_starts)) + ','])
        return s
    
#    def to_gtf(self):
#        f = GTFFeature()
#        f.seqid = self.chrom
#        f.source = 'assemblyline'
#        f.feature_type = 'transcript'
#        f.start = self.start
#        f.end = self.end
#        f.score = 1000.0
#        f.strand = strand_int_to_str(self.strand)
#        f.phase = '.'
#        f.attrs = {}
#        self.id = None
#        self.library = None
#        self.sample = None
#        self.cohort = None
#        self.fpkm = 0.0
#        self.cov = 0.0
#        self.frac = 0.0
#        self.exons = None
#        self.attrs = None

    
    
    
