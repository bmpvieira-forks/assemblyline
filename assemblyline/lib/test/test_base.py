'''
Created on Feb 20, 2011

@author: mkiyer
'''
import os
import networkx as nx

from ..transcript_graph import TranscriptGraph
from ..transcript_parser import Transcript, parse_gtf
from ..base import Exon, strand_str_to_int

DOT_DIR = "dot_files"
GTF_DIR = "gtf_files"

def make_transcript(exons,
                    chrom='chr1',
                    id='G1.1',
                    strand='.',
                    score=1):
    start = min(e[0] for e in exons)
    end = max(e[1] for e in exons)
    length = sum(e[1] - e[0] for e in exons)
    t = Transcript()
    t.id = id
    t.label = None
    t.strand = strand_str_to_int(strand)
    t.score = float(score)
    t.chrom = chrom
    t.start = start
    t.end = end
    t.length = length
    t.exons = [Exon(e[0], e[1]) for e in exons]
    return t

def get_dot_path(filename):
    return os.path.join(os.path.dirname(__file__), DOT_DIR, filename)

def get_gtf_path(filename):
    return os.path.join(os.path.dirname(__file__), GTF_DIR, filename)

def convert_data_to_string(exondata):
    return ','.join([str(d) for d in exondata])    

def convert_attrs_to_strings(H):
    G = H.copy()
    for n in G:
        G.node[n]['data'] = convert_data_to_string(G.node[n]['data'])
    for u,v,d in G.edges_iter(data=True):
        d['data'] = convert_data_to_string(d['data'])
    return G

def compare_dot(transcript_graph, filename):
    Gcorrect = nx.read_dot(get_dot_path(filename))
    G = convert_attrs_to_strings(transcript_graph.G)
    for n in G:
        node_str = str(n)
        if node_str not in Gcorrect:
            return False
        correct_data = Gcorrect.node[node_str]['data']
        test_data = '"%s"' % (G.node[n]['data'])
        if correct_data != test_data:
            return False
    for u,v,d in G.edges_iter(data=True):
        ustr = str(u)
        vstr = str(v)
        correct_data = Gcorrect.edge[ustr][vstr][0]['data']
        test_data = '"%s"' % d['data']
        if correct_data != test_data:
            return False
    return nx.is_isomorphic(transcript_graph.G, Gcorrect)

def write_dot(transcript_graph, filename):
    G = convert_attrs_to_strings(transcript_graph.G)
    nx.write_dot(G, get_dot_path(filename))

def read_gtf(filename):
    for locus_transcripts in parse_gtf(open(get_gtf_path(filename))):
        isoform_graph = TranscriptGraph.from_transcripts(locus_transcripts)
        return isoform_graph