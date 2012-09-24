'''
Created on Dec 17, 2011

@author: mkiyer
'''
from assemblyline.lib.transcript import Exon

# transcript graph attributes
STRAND_SCORE = 'strand_score'
TRANSCRIPT_IDS = 'ids'

# assembler graph attributes
NODE_SCORE = 'score'
NODE_LENGTH = 'length'
SMOOTH_FWD = 'smfwd'
SMOOTH_REV = 'smrev'
SMOOTH_TMP = 'smtmp'
CHAIN_NODES = 'chain'
CHAIN_DATA = 'chain_data'
CHAIN_EDGES = 'chain_edges'

# constant minimum path score
MIN_SCORE = 1.0e-8

# global variables for current locus, gene, tss, transcript id
GLOBAL_LOCUS_ID = 1
GLOBAL_GENE_ID = 1
GLOBAL_TSS_ID = 1
GLOBAL_TRANSCRIPT_ID = 1

def init_directed_node_attrs(n):
    """
    each node in directed graph maintains attributes:
    'ids': set() of transcript id strings
    'score': node weight (float)
    'length': size of node in nucleotides
    """
    return {TRANSCRIPT_IDS: set(),
            NODE_LENGTH: (n.end - n.start), 
            NODE_SCORE: 0.0}
