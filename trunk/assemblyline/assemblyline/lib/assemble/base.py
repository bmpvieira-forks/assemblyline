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

# constant minimum path density
MIN_DENSITY = 1.0e-8

# global variables for current locus, gene, tss, transcript id
GLOBAL_LOCUS_ID = 1
GLOBAL_GENE_ID = 1
GLOBAL_TSS_ID = 1
GLOBAL_TRANSCRIPT_ID = 1