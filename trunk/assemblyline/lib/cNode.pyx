'''
Created on Nov 29, 2010

@author: mkiyer
'''
cimport cython

DEF POS_STRAND = 0
DEF NEG_STRAND = 1
DEF NO_STRAND = 2
cdef list strand_int_to_str = ["+", "-", "."]

DEF EXON = 0
DEF INTRON = 1
cdef list node_type_to_str = ["E", "I"]

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
        return ("%s-%s(%s)[%s]" % (self.start, self.end, strand_int_to_str[self.strand], node_type_to_str[self.node_type][0]))
    def __repr__(Node self):
        return ("<%s(start=%d end=%d strand=%s node_type=%s>" %
                (self.__class__.__name__, self.start, self.end, strand_int_to_str[self.strand], node_type_to_str[self.node_type]))

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
        if op == 2:
            return eq
        elif op == 3:
            return not eq
        else:
            print "op is", op
            assert False
        
    def __hash__(self):
        return (self.start << 18) | (self.end << 3) | (self.strand << 1) | (self.node_type)

    def cmp_strand(self, other):
        if self.strand == other.strand:
            return True
        elif (self.strand == NO_STRAND) or (other.strand == NO_STRAND):
            return True
        return False

    def is_overlapping(self, other):
        return (self.start < other.end) and (other.start < self.end) 
