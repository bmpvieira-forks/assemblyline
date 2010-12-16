'''
Created on Nov 29, 2010

@author: mkiyer
'''
cimport cython

DEF POS_STRAND = 0
DEF NEG_STRAND = 1
DEF NO_STRAND = 2

cdef list _strand_int_to_str = ["+", "-", "."]
cdef dict _strand_str_to_int = {"+": POS_STRAND,
                                "-": NEG_STRAND,
                                ".": NO_STRAND}

DEF EXON = 0
DEF INTRON = 1
DEF DUMMY = 2

cdef list _node_type_to_str = ["E", "I", "D"]

def strand_str_to_int(strand):
    return _strand_str_to_int[strand]
def strand_int_to_str(strand):
    return _strand_int_to_str[strand]
def node_type_to_str(int node_type):
    return _node_type_to_str[node_type]

cdef class Node(object):
    
    cdef public int start
    cdef public int end
    cdef public int strand
    cdef public int node_type

    def __cinit__(Node self):
        self.start = 0
        self.end = 0
        self.strand = NO_STRAND 
        self.node_type = EXON
    
    def __init__(Node self, int start, int end, int strand, int node_type):
        self.start = start
        self.end = end
        self.strand = strand
        self.node_type = node_type
    
    def __deepcopy__(Node other, object memo):        
        return Node(other.start, other.end, other.strand, other.node_type)
    
    def __str__(Node self):
        return ("%s-%s(%s)[%s]" % (self.start, self.end, _strand_int_to_str[self.strand], _node_type_to_str[self.node_type][0]))
    def __repr__(Node self):
        return ("<%s(start=%d end=%d strand=%s node_type=%s>" %
                (self.__class__.__name__, self.start, self.end, _strand_int_to_str[self.strand], _node_type_to_str[self.node_type]))

    def __richcmp__(Node self, Node other, int op):
        """
        <     0
        ==     2
        >     4
        <=     1
        !=     3
        >=     5
        """
        cdef bint eq
        eq = ((self.start == other.start) and
              (self.end == other.end) and
              (self.strand == other.strand) and
              (self.node_type == other.node_type))
        if op == 0:
            return self.start < other.start 
        elif op == 2:
            return eq
        elif op == 3:
            return not eq
        else:
            print "op is", op
            assert False
        
    def __hash__(Node self):
        return (self.start << 18) | (self.end << 3) | (self.strand << 1) | (self.node_type)

    def cmp_strand(Node self, Node other):
        return cmp_strand(self.strand, other.strand)

    def is_overlapping(self, other):
        return interval_overlap(self, other)

cpdef bint cmp_strand(int a, int b):
    if (a == NO_STRAND) or (b == NO_STRAND):
        return True
    return a == b

cpdef bint interval_overlap(Node a, Node b):
    return (a.start < b.end) and (b.start < a.end)

cpdef bint interval_overlap_threshold(Node a, Node b, int d=0):
    return (a.start < b.end + d) and (b.start < a.end + d)
    
