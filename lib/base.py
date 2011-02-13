'''
Created on Dec 4, 2010

@author: mkiyer
'''
EXON = 0
INTRON = 1
DUMMY = 2
POS_STRAND = 0
NEG_STRAND = 1
NO_STRAND = 2

def merge_strand(strand1, strand2):
    # TODO: move to C code
    if strand1 == strand2:
        return strand1
    elif strand1 == NO_STRAND:
        return strand2
    elif strand2 == NO_STRAND:
        return strand1
    else:
        assert False