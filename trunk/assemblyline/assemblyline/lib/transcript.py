'''
Created on Nov 30, 2011

@author: mkiyer
'''
from gtf import GTFFeature

# attributes
TRANSCRIPT_ID = "transcript_id"
GENE_ID = "gene_id"

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

def merge_strand(strand1, strand2):
    if strand1 == strand2:
        return strand1
    elif strand1 == NO_STRAND:
        return strand2
    elif strand2 == NO_STRAND:
        return strand1
    else:
        assert False

def cmp_strand(a, b):
    if (a == NO_STRAND) or (b == NO_STRAND):
        return True
    return a == b

def interval_overlap(a, b):
    return (a.start < b.end) and (b.start < a.end)

def interval_overlap_threshold(a, b, d=0):
    return (a.start < b.end + d) and (b.start < a.end + d)

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
    __slots__ = ('chrom', 'start', 'end', 'strand', 'score', 'exons', 'attrs') 

    def __init__(self):
        self.chrom = None
        self.start = -1
        self.end = -1
        self.strand = NO_STRAND
        self.score = 0.0
        self.exons = None
        self.attrs = {}

    def __str__(self):
        return ("<%s(chrom='%s', start='%d', end='%d', strand='%s', "
                "score='%s' exons='%s', attrs='%s'" %
                (self.__class__.__name__, self.chrom, self.start, self.end, 
                 strand_int_to_str(self.strand), str(self.score), self.exons, 
                 self.attrs))
    @property
    def length(self):
        return sum((e.end - e.start) for e in self.exons)

    def iterintrons(self):
        e1 = self.exons[0]
        for e2 in self.exons[1:]:
            yield (e1.end,e2.start)
            e1 = e2

    def to_bed12(self):
        block_sizes = []
        block_starts = []
        for e0, e1 in self.exons:
            block_starts.append(e0 - self.tx_start)
            block_sizes.append(e1 - e0)        
        # write
        s = '\t'.join([self.chrom, 
                       str(self.start), 
                       str(self.end),
                       str(self.attrs["transcript_id"]),
                       '0',
                       strand_int_to_str(self.strand),
                       str(self.start),
                       str(self.start),
                       '0',
                       str(len(self.exons)),
                       ','.join(map(str,block_sizes)) + ',',
                       ','.join(map(str,block_starts)) + ','])
        return s
    
    def to_gtf_features(self):
        # transcript feature
        f = GTFFeature()
        f.seqid = self.chrom
        f.source = 'assemblyline'
        f.feature_type = 'transcript'
        f.start = self.start
        f.end = self.end
        f.score = 1000.0
        f.strand = strand_int_to_str(self.strand)
        f.phase = '.'
        f.attrs = self.attrs
        features = [f]
        # exon features
        for i,e in enumerate(self.exons):
            f = GTFFeature()
            f.seqid = self.chrom
            f.source = 'assemblyline'
            f.feature_type = 'exon'
            f.start = e.start
            f.end = e.end
            f.score = 1000.0
            f.strand = strand_int_to_str(self.strand)
            f.phase = '.'
            f.attrs = {}
            f.attrs["gene_id"] = self.attrs["gene_id"]
            f.attrs["transcript_id"] = self.attrs["transcript_id"]
            f.attrs["exon_number"] = i
            features.append(f)
        return features
