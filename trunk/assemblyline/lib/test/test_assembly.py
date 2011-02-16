'''
Created on Feb 16, 2011

@author: mkiyer
'''
import logging
import unittest
import os
import networkx as nx
import matplotlib.pyplot as plt

from ..isoform_graph import IsoformGraph, get_in_degree, get_out_degree
from ..transcript_parser import Transcript, Exon
from ..cnode import Node, strand_str_to_int
from ..base import POS_STRAND, NEG_STRAND, NO_STRAND

logging.basicConfig(level=logging.DEBUG,
                    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")


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
    t.score = score
    t.chrom = chrom
    t.start = start
    t.end = end
    t.length = length
    t.exons = [Exon(e[0], e[1]) for e in exons]
    return t

class Test(unittest.TestCase):

    def test_strand_degree(self):
        """ensure that in_degree/out_degree are correctly ignoring nodes
        on the opposite strand"""        
        ig = IsoformGraph()
        t1 = make_transcript(((0, 100), (500,600), (900, 1000)), strand="+")
        t2 = make_transcript(((300, 400), (500,600), (700,800)), strand="-")
        ig.add_transcripts([t1, t2])
        G = ig.G        
        self.assertTrue(get_in_degree(G, n, POS_STRAND)

        for n in G.nodes_iter():
            print n, get_in_degree(G, n, 0), get_out_degree(G, n, 0)




if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()