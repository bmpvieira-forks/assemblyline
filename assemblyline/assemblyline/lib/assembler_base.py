'''
Created on Dec 3, 2011

@author: mkiyer
'''
from transcript import Exon

# assembler graph attributes
NODE_DENSITY = 'density'
EDGE_OUT_FRAC = 'outfrac'
EDGE_IN_FRAC = 'infrac'
EDGE_DENSITY = 'density'

# constants for fake 'start' and 'end' nodes
SOURCE_NODE = Exon(-1,-1)
SINK_NODE = Exon(-2,-2)
