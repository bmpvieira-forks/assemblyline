'''
Created on Feb 14, 2011

@author: mkiyer
'''
import logging
import unittest
import os
import networkx as nx
import matplotlib.pyplot as plt
from StringIO import StringIO

from ..base import Exon, POS_STRAND, NEG_STRAND
from ..transcript_graph import TranscriptGraph
from ..transcript_parser import parse_gtf

from test_base import compare_dot, write_dot, convert_attrs_to_strings, read_gtf, make_transcript

logging.basicConfig(level=logging.DEBUG,
                    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")

class TestTranscriptGraph(unittest.TestCase):

    def test_unstranded_transcripts(self):
        """
        ensure that coverage from unstranded transcripts is allocated
        proportionally to the coverage from overlapping stranded 
        transcripts
        """
        tg = TranscriptGraph()
        t1 = make_transcript(((0, 100), (500,600), (900, 1100)), id="1", strand="+", score=20)
        t2 = make_transcript(((0, 100), (500,600), (900, 1100)), id="1", strand="-", score=20)
        t3 = make_transcript([(1050, 1150), (1500,1600)], id="2", strand=".", score=20)
        # ensure no negative strand is part of graph
        tg.add_transcripts([t1, t3])
        for n,d in tg.G.nodes_iter(data=True):
            self.assertTrue(all(x.scores[NEG_STRAND] == 0 for x in d['data']))
        # ensure no positive strand is part of graph except for the
        # final exon that does not overlap with a negative strand exon
        tg.add_transcripts([t2, t3])
        for u,v,d in tg.G.edges_iter(data=True):
            self.assertTrue(d['strand'] != POS_STRAND)
        self.assertTrue(tg.G.node[Exon(1500,1600)]['data'][0].scores == [10.0, 0.0])
        # now combine both plus and minus strands
        t1 = make_transcript(((0, 100), (500,600)), id="1", strand="+", score=20)
        t2 = make_transcript(((500,900),), id="2", strand=".", score=20)
        t3 = make_transcript([(700,900), (1500,1700)], id="3", strand="-", score=20)
        # ensure coverage is partitioned to + and - strand
        tg.add_transcripts([t1, t2, t3])
        self.assertTrue(tg.G.node[Exon(600,700)]['data'][0].scores == [2.5,2.5])        
        #nx.spring_layout(G)
        #nx.draw(G)
        #plt.show()

    def test_basic_graph(self):
        """test that basic transcript parser works and builds graph"""        
        # plus strand
        test_basename = "basic_graph1"
        gtf_file = test_basename + ".gtf"
        dot_file = test_basename + ".dot"
        txgraph = read_gtf(gtf_file)
        #write_dot(txgraph, dot_file)
        self.assertTrue(compare_dot(txgraph, dot_file))
        # minus strand
        test_basename = "basic_graph2"
        gtf_file = test_basename + ".gtf"
        dot_file = test_basename + ".dot"
        txgraph = read_gtf(gtf_file)
        #write_dot(txgraph, dot_file)
        self.assertTrue(compare_dot(txgraph, dot_file))
        # unstranded multi-exon
        test_basename = "basic_graph3"
        gtf_file = test_basename + ".gtf"
        dot_file = test_basename + ".dot"
        txgraph = read_gtf(gtf_file)
        #write_dot(txgraph, dot_file)
        self.assertTrue(compare_dot(txgraph, dot_file))
        #nx.spring_layout(txgraph.G)
        #nx.draw(txgraph.G)
        #plt.show()

    def test_basic_split(self):
        """test that initial graph construction correctly splits exons"""
        # plus strand
        test_basename = "basic_split1"
        gtf_file = test_basename + ".gtf"
        dot_file = test_basename + ".dot"
        txgraph = read_gtf(gtf_file)
        #write_dot(txgraph, dot_file)
        compare_dot(txgraph, dot_file)
        # minus strand
        test_basename = "basic_split2"
        gtf_file = test_basename + ".gtf"
        dot_file = test_basename + ".dot"
        txgraph = read_gtf(gtf_file)
        #write_dot(txgraph, dot_file)
        compare_dot(txgraph, dot_file)
        # unstranded
        test_basename = "basic_split3"
        gtf_file = test_basename + ".gtf"
        dot_file = test_basename + ".dot"
        txgraph = read_gtf(gtf_file)
        #write_dot(txgraph, dot_file)
        compare_dot(txgraph, dot_file)
        #nx.spring_layout(txgraph.G)
        #nx.draw(txgraph.G)
        #plt.show()

    def test_complex_split(self):
        """test some more complicated exon splitting cases"""
        # only one strand
        test_basename = "complex_split1"
        gtf_file = test_basename + ".gtf"
        dot_file = test_basename + ".dot"
        txgraph = read_gtf(gtf_file)
        #write_dot(txgraph, dot_file)        
        compare_dot(txgraph, dot_file)
        # duplicate plus and minus strand overlapping
        test_basename = "complex_split2"
        gtf_file = test_basename + ".gtf"
        dot_file = test_basename + ".dot"
        txgraph = read_gtf(gtf_file)
        #write_dot(txgraph, dot_file)
        compare_dot(txgraph, dot_file)
        # add in unstranded transcripts
        test_basename = "complex_split3"
        gtf_file = test_basename + ".gtf"
        dot_file = test_basename + ".dot"
        txgraph = read_gtf(gtf_file)
        #write_dot(txgraph, dot_file)
        compare_dot(txgraph, dot_file)
        #nx.spring_layout(txgraph.G)
        #nx.draw(txgraph.G)
        #plt.show()

    def test_antisense(self):
        """test some examples of antisense transcripts"""
        # exons do not overlap
        test_basename = "antisense1"
        gtf_file = test_basename + ".gtf"
        dot_file = test_basename + ".dot"
        txgraph = read_gtf(gtf_file)
        #write_dot(txgraph, dot_file)
        compare_dot(txgraph, dot_file)
        # exons partially overlap
        test_basename = "antisense2"
        gtf_file = test_basename + ".gtf"
        dot_file = test_basename + ".dot"
        txgraph = read_gtf(gtf_file)
        #write_dot(txgraph, dot_file)
        compare_dot(txgraph, dot_file)
        # one exon shared among strands
        test_basename = "antisense3"
        gtf_file = test_basename + ".gtf"
        dot_file = test_basename + ".dot"
        txgraph = read_gtf(gtf_file)
        #write_dot(txgraph, dot_file)
        compare_dot(txgraph, dot_file)




if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()