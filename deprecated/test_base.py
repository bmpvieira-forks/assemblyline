'''
Created on Feb 20, 2011

@author: mkiyer
'''
import os
import networkx as nx

from assemblyline.lib.transcript import Exon, strand_str_to_int
from assemblyline.lib import transcript_graph
from assemblyline.lib.transcript_parser import Transcript, parse_gtf, cufflinks_attr_defs

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

def convert_attrs_to_strings(H):
    G = H.copy()
    for n,d in G.nodes_iter(data=True):
        d['ids'] = str(sorted(d['ids']))
        d['strand_density'] = str(['%.3f' % x for x in d['strand_density']])
    return G

def compare_dot(H, filename):
    Gcorrect = nx.read_dot(get_dot_path(filename))
    G = convert_attrs_to_strings(H)
    for n,d in G.nodes_iter(data=True):
        node_str = str(n)
        if node_str not in Gcorrect:
            return False
        correct_ids = Gcorrect.node[node_str]['ids'].strip('"')
        if correct_ids != d['ids']:
            return False
        correct_strand_density = Gcorrect.node[node_str]['strand_density'].strip('"')
        if correct_strand_density != d['strand_density']:
            return False
    for u,v,d in G.edges_iter(data=True):
        ustr = str(u)
        vstr = str(v)
        if not Gcorrect.has_edge(ustr, vstr):
            return False
    return nx.is_isomorphic(H, Gcorrect)

def write_dot(H, filename):
    G = convert_attrs_to_strings(H)
    nx.write_dot(G, get_dot_path(filename))

def read_gtf(filename):
    """
    reads and returns graph from first locus in GTF file
    """
    for locus_transcripts in parse_gtf(open(get_gtf_path(filename)), cufflinks_attr_defs):
        G = transcript_graph.add_transcripts_to_graph(locus_transcripts)
        return G
