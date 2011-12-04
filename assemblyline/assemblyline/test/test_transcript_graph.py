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

from assemblyline.lib.transcript import Exon, POS_STRAND, NEG_STRAND
from assemblyline.lib.transcript_graph import trim_left, trim_right, trim_transcript
from assemblyline.lib.transcript_parser import parse_gtf

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
        #nx.spring_layout(txgraph.G)
        #nx.draw(txgraph.G)
        #plt.show()
       
#    def test_example1(self):
#        """test some examples from actual data"""
#        # exons do not overlap
#        test_basename = "ex1"
#        gtf_file = test_basename + ".gtf"
#        dot_file = test_basename + ".dot"
#        txgraph = read_gtf(gtf_file)
#        write_dot(txgraph, dot_file)
#        self.assertTrue(compare_dot(txgraph, dot_file))        


class TestTrimming(unittest.TestCase):

    def test_trim_left_right(self):
        # make up some intron start/end intervals
        intron_starts = [100, 400, 700]
        intron_ends = [200, 500, 800]
        # try trimming fake exon left
        overhang_threshold = 10
        x = trim_left(195, 220, intron_starts, intron_ends, overhang_threshold)
        self.assertEqual(x, 200)
        x = trim_left(191, 220, intron_starts, intron_ends, overhang_threshold)
        self.assertEqual(x, 200)
        x = trim_left(190, 220, intron_starts, intron_ends, overhang_threshold)
        self.assertEqual(x, 200)
        x = trim_left(189, 220, intron_starts, intron_ends, overhang_threshold)
        self.assertEqual(x, 189)
        # try trimming fake exon right
        x = trim_right(300, 405, intron_starts, intron_ends, overhang_threshold)
        self.assertEqual(x, 400)
        x = trim_right(300, 409, intron_starts, intron_ends, overhang_threshold)
        self.assertEqual(x, 400)
        x = trim_right(300, 410, intron_starts, intron_ends, overhang_threshold)
        self.assertEqual(x, 400)
        x = trim_right(300, 411, intron_starts, intron_ends, overhang_threshold)
        self.assertEqual(x, 411)
        
    def test_trim_transcript(self):        
        # make up some intron start/end intervals
        intron_starts = [100, 400, 700]
        intron_ends = [200, 500, 800]
        overhang_threshold = 10        
        # trim some fake transcripts
        s,e = trim_transcript([Exon(195,405)], intron_starts, intron_ends, 
                              overhang_threshold)
        self.assertEquals((s,e), (200,400))
        s,e = trim_transcript([Exon(190,410)], intron_starts, intron_ends, 
                              overhang_threshold)
        self.assertEquals((s,e), (200,400))
        s,e = trim_transcript([Exon(189,411)], intron_starts, intron_ends, 
                              overhang_threshold)
        self.assertEquals((s,e), (189,411))
        # create situation where transcript trimming conflicts
        intron_starts = [100, 190, 400]
        intron_ends = [200, 500]
        overhang_threshold = 20
        s,e = trim_transcript([Exon(188,201)], intron_starts, intron_ends, 
                              overhang_threshold)
        self.assertEquals((s,e), (188,190))
        # behavior in the case of a tie is to trim right (arbitrary)
        s,e = trim_transcript([Exon(189,201)], intron_starts, intron_ends, 
                              overhang_threshold)
        self.assertEquals((s,e), (200,201))
        s,e = trim_transcript([Exon(189,202)], intron_starts, intron_ends, 
                              overhang_threshold)
        self.assertEquals((s,e), (200,202))

#    def test_trimming(self):
#        """test some examples of trimming transcripts"""
#        # simple example
#        test_basename = "trim1"
#        gtf_file = test_basename + ".gtf"
#        dot_file = test_basename + ".dot"
#        # without trimming
#        txgraph = read_gtf(gtf_file)
#        #write_dot(txgraph, dot_file)
#        self.assertTrue(compare_dot(txgraph, dot_file))        
#        # with trimming
#        dot_file = test_basename + "_trim15.dot"
#        for locus_transcripts in parse_gtf(open(get_gtf_path(gtf_file))):
#            txgraph = TranscriptGraph()
#            txgraph.add_transcripts(locus_transcripts, overhang_threshold=15)
#            break
#        #write_dot(txgraph, dot_file)
#        self.assertTrue(compare_dot(txgraph, dot_file))        


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()