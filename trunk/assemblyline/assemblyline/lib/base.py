'''
Created on Feb 15, 2011

@author: mkiyer
'''

# strand constants
POS_STRAND = 0
NEG_STRAND = 1
NO_STRAND = 2

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

