'''
Created on Feb 14, 2011

@author: mkiyer
'''
import unittest
import os
import networkx as nx
import matplotlib.pyplot as plt

from ..isoform_graph import IsoformGraph
from ..transcript_parser import parse_gtf

import logging
logging.basicConfig(level=logging.DEBUG,
                    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")

DOT_DIR = "dot_files"
GTF_DIR = "gtf_files"

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
    return G

def write_dot(isoform_graph, filename):
    G = convert_attrs_to_strings(isoform_graph.G)
    nx.write_dot(G, get_dot_path(filename))

def read_gtf(filename):
    for locus_transcripts in parse_gtf(open(get_gtf_path(filename))):
        print locus_transcripts
        isoform_graph = IsoformGraph.from_transcripts(locus_transcripts)
        return isoform_graph

class TestIsoformGraph(unittest.TestCase):

    def compare_dot(self, isoform_graph, filename):
        Gcorrect = nx.read_dot(get_dot_path(filename))
        G = convert_attrs_to_strings(isoform_graph.G)
        for n in G:
            node_str = str(n)
            self.assertTrue(node_str in Gcorrect)
            correct_data = Gcorrect.node[node_str]['data']
            test_data = '"%s"' % (G.node[n]['data'])
            self.assertTrue(correct_data == test_data)
        self.assertTrue(nx.is_isomorphic(isoform_graph.G, Gcorrect)) 

    def testBasicGraph(self):
        """test that basic transcript parser works and builds graph"""        
        # plus strand
        test_basename = "basic_graph1"
        gtf_file = test_basename + ".gtf"
        dot_file = test_basename + ".dot"
        isoform_graph = read_gtf(gtf_file)
        #write_dot(isoform_graph, dot_file)
        self.compare_dot(isoform_graph, dot_file)
        # minus strand
        test_basename = "basic_graph2"
        gtf_file = test_basename + ".gtf"
        dot_file = test_basename + ".dot"
        isoform_graph = read_gtf(gtf_file)
        #write_dot(isoform_graph, dot_file)
        self.compare_dot(isoform_graph, dot_file)
        # unstranded multi-exon
        test_basename = "basic_graph3"
        gtf_file = test_basename + ".gtf"
        dot_file = test_basename + ".dot"
        isoform_graph = read_gtf(gtf_file)
        #write_dot(isoform_graph, dot_file)
        self.compare_dot(isoform_graph, dot_file)
        #nx.spring_layout(isoform_graph.G)
        #nx.draw(isoform_graph.G)
        #plt.show()

    def testBasicSplit(self):
        """test that initial graph construction correctly splits exons"""
        # plus strand
        test_basename = "basic_split1"
        gtf_file = test_basename + ".gtf"
        dot_file = test_basename + ".dot"
        isoform_graph = read_gtf(gtf_file)
        #write_dot(isoform_graph, dot_file)
        self.compare_dot(isoform_graph, dot_file)
        # minus strand
        test_basename = "basic_split2"
        gtf_file = test_basename + ".gtf"
        dot_file = test_basename + ".dot"
        isoform_graph = read_gtf(gtf_file)
        #write_dot(isoform_graph, dot_file)
        self.compare_dot(isoform_graph, dot_file)
        # unstranded
        test_basename = "basic_split3"
        gtf_file = test_basename + ".gtf"
        dot_file = test_basename + ".dot"
        isoform_graph = read_gtf(gtf_file)
        #write_dot(isoform_graph, dot_file)
        self.compare_dot(isoform_graph, dot_file)
        #nx.spring_layout(isoform_graph.G)
        #nx.draw(isoform_graph.G)
        #plt.show()

    def testComplexSplit(self):
        """test some more complicated exon splitting cases"""
        # only one strand
        test_basename = "complex_split1"
        gtf_file = test_basename + ".gtf"
        dot_file = test_basename + ".dot"
        isoform_graph = read_gtf(gtf_file)
        #write_dot(isoform_graph, dot_file)        
        self.compare_dot(isoform_graph, dot_file)
        # duplicate plus and minus strand overlapping
        test_basename = "complex_split2"
        gtf_file = test_basename + ".gtf"
        dot_file = test_basename + ".dot"
        isoform_graph = read_gtf(gtf_file)
        #write_dot(isoform_graph, dot_file)
        self.compare_dot(isoform_graph, dot_file)
        # add in unstranded transcripts
        test_basename = "complex_split3"
        gtf_file = test_basename + ".gtf"
        dot_file = test_basename + ".dot"
        isoform_graph = read_gtf(gtf_file)
        #write_dot(isoform_graph, dot_file)
        self.compare_dot(isoform_graph, dot_file)
        #nx.spring_layout(isoform_graph.G)
        #nx.draw(isoform_graph.G)
        #plt.show()

    def testAntisense(self):
        """test some examples of antisense transcripts"""
        # exons do not overlap
        test_basename = "antisense1"
        gtf_file = test_basename + ".gtf"
        dot_file = test_basename + ".dot"
        isoform_graph = read_gtf(gtf_file)
        #write_dot(isoform_graph, dot_file)
        self.compare_dot(isoform_graph, dot_file)
        # exons partially overlap
        test_basename = "antisense2"
        gtf_file = test_basename + ".gtf"
        dot_file = test_basename + ".dot"
        isoform_graph = read_gtf(gtf_file)
        #write_dot(isoform_graph, dot_file)
        self.compare_dot(isoform_graph, dot_file)
        # one exon shared among strands
        test_basename = "antisense3"
        gtf_file = test_basename + ".gtf"
        dot_file = test_basename + ".dot"
        isoform_graph = read_gtf(gtf_file)
        #write_dot(isoform_graph, dot_file)
        self.compare_dot(isoform_graph, dot_file)


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()