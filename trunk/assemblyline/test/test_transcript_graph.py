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

from lib.base import Exon, POS_STRAND, NEG_STRAND
from lib.transcript_graph import TranscriptGraph
from lib.transcript_parser import parse_gtf

from test_base import compare_dot, write_dot, convert_attrs_to_strings, read_gtf, make_transcript, get_gtf_path

logging.basicConfig(level=logging.DEBUG,
                    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")

class TestTranscriptGraph(unittest.TestCase):

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
        self.assertTrue(compare_dot(txgraph, dot_file))
        # minus strand
        test_basename = "basic_split2"
        gtf_file = test_basename + ".gtf"
        dot_file = test_basename + ".dot"
        txgraph = read_gtf(gtf_file)
        #write_dot(txgraph, dot_file)
        self.assertTrue(compare_dot(txgraph, dot_file))
        # unstranded
        test_basename = "basic_split3"
        gtf_file = test_basename + ".gtf"
        dot_file = test_basename + ".dot"
        txgraph = read_gtf(gtf_file)
        #write_dot(txgraph, dot_file)
        self.assertTrue(compare_dot(txgraph, dot_file))
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
        self.assertTrue(compare_dot(txgraph, dot_file))
        # duplicate plus and minus strand overlapping
        test_basename = "complex_split2"
        gtf_file = test_basename + ".gtf"
        dot_file = test_basename + ".dot"
        txgraph = read_gtf(gtf_file)
        #write_dot(txgraph, dot_file)
        self.assertTrue(compare_dot(txgraph, dot_file))
        # add in unstranded transcripts
        test_basename = "complex_split3"
        gtf_file = test_basename + ".gtf"
        dot_file = test_basename + ".dot"
        txgraph = read_gtf(gtf_file)
        #write_dot(txgraph, dot_file)
        self.assertTrue(compare_dot(txgraph, dot_file))
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
        self.assertTrue(compare_dot(txgraph, dot_file))        
        # exons partially overlap
        test_basename = "antisense2"
        gtf_file = test_basename + ".gtf"
        dot_file = test_basename + ".dot"
        txgraph = read_gtf(gtf_file)
        #write_dot(txgraph, dot_file)
        self.assertTrue(compare_dot(txgraph, dot_file))        
        # one exon shared among strands
        test_basename = "antisense3"
        gtf_file = test_basename + ".gtf"
        dot_file = test_basename + ".dot"
        txgraph = read_gtf(gtf_file)
        #write_dot(txgraph, dot_file)
        self.assertTrue(compare_dot(txgraph, dot_file))        

#    def test_example1(self):
#        """test some examples from actual data"""
#        pass
#        # exons do not overlap
#        #test_basename = "ex1"
#        #gtf_file = test_basename + ".gtf"
#        #dot_file = test_basename + ".dot"
#        #txgraph = read_gtf(gtf_file)
#        #write_dot(txgraph, dot_file)
#        #self.assertTrue(compare_dot(txgraph, dot_file))        


class TestTrimming(unittest.TestCase):

    def test_trimming(self):
        """test some examples of trimming transcripts"""
        # simple example
        test_basename = "trim1"
        gtf_file = test_basename + ".gtf"
        dot_file = test_basename + ".dot"
        # without trimming
        txgraph = read_gtf(gtf_file)
        #write_dot(txgraph, dot_file)
        self.assertTrue(compare_dot(txgraph, dot_file))        
        # with trimming
        dot_file = test_basename + "_trim15.dot"
        for locus_transcripts in parse_gtf(open(get_gtf_path(gtf_file))):
            txgraph = TranscriptGraph()
            txgraph.add_transcripts(locus_transcripts, overhang_threshold=15)
            break
        #write_dot(txgraph, dot_file)
        self.assertTrue(compare_dot(txgraph, dot_file))        


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()