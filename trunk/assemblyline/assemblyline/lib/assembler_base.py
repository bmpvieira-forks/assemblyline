'''
Created on Dec 3, 2011

@author: mkiyer
'''
from transcript import Exon

# transcript graph attributes
STRAND_DENSITY = 'strand_density'
TRANSCRIPT_IDS = 'ids'

# assembler graph attributes
NODE_DENSITY = 'density'
NODE_LENGTH = 'length'
EDGE_OUT_FRAC = 'outfrac'
EDGE_IN_FRAC = 'infrac'

# constants for fake 'start' and 'end' nodes
SOURCE_NODE = Exon(-1,-1)
SINK_NODE = Exon(-2,-2)

# constant minimum node density
MIN_NODE_DENSITY = 1.0e-6

# global variables for current locus, gene, tss, transcript id
GLOBAL_LOCUS_ID = 1
GLOBAL_GENE_ID = 1
GLOBAL_TSS_ID = 1
GLOBAL_TRANSCRIPT_ID = 1
