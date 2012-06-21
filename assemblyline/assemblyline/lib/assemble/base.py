'''
Created on Dec 17, 2011

@author: mkiyer
'''
from assemblyline.lib.transcript import Exon

# transcript graph attributes
STRAND_DENSITY = 'strand_density'
TRANSCRIPT_IDS = 'ids'

# assembler graph attributes
NODE_DENSITY = 'density'
NODE_LENGTH = 'length'
NODE_TSS_ID = 'tss_id'

# constants for fake 'start' and 'end' nodes
SOURCE_NODE = Exon(-1,-1)
SINK_NODE = Exon(-2,-2)

# constants for direction in path growing
FWD = 0
REV = 1

# constant minimum path density
MIN_DENSITY = 1.0e-8

# global variables for current locus, gene, tss, transcript id
GLOBAL_LOCUS_ID = 1
GLOBAL_GENE_ID = 1
GLOBAL_TSS_ID = 1
GLOBAL_TRANSCRIPT_ID = 1