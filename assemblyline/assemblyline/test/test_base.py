'''
Created on Nov 8, 2012

@author: mkiyer
'''
import os

from assemblyline.lib.transcript import parse_gtf
from assemblyline.lib.assemble.transcript_graph import create_transcript_graphs
from assemblyline.lib.base import GTFAttr

GTF_DIR = "gtf_files"

def get_gtf_path(filename):
    return os.path.join(os.path.dirname(__file__), GTF_DIR, filename)

def read_gtf(filename):
    return list(parse_gtf(open(get_gtf_path(filename))))

def read_first_locus(filename):
    loci = read_gtf(filename)
    for t in loci[0]:
        t.score = float(t.attrs[GTFAttr.SCORE])
    return loci[0]

def get_transcript_graphs(transcripts):
    GG = {}
    for G, strand, strand_transcript_map in \
        create_transcript_graphs(transcripts):
        GG[strand] = (G, strand_transcript_map)
    return GG