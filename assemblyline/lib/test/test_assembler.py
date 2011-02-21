'''
Created on Feb 16, 2011

@author: mkiyer
'''
import logging
import unittest
import os
import networkx as nx
import matplotlib.pyplot as plt

from ..base import POS_STRAND, NEG_STRAND, NO_STRAND
from ..transcript_graph import TranscriptGraph, ExonData, get_in_degree, get_out_degree
from ..assembler import calculate_node_weight

from test_base import make_transcript, read_gtf, write_dot

logging.basicConfig(level=logging.DEBUG,
                    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")

class TestAssembler(unittest.TestCase):

    def test_calculate_node_weight(self):
        # test all positive strand nodes
        exon_data_list = [ExonData("1", POS_STRAND, 10),
                          ExonData("2", POS_STRAND, 20)]
        weights = calculate_node_weight(exon_data_list)
        self.assertTrue(weights == (30,0))
        # test all negative strand nodes        
        exon_data_list = [ExonData("1", NEG_STRAND, 5),
                          ExonData("2", NEG_STRAND, 6)]
        weights = calculate_node_weight(exon_data_list)
        self.assertTrue(weights == (0,11))
        # test only unstranded nodes
        exon_data_list = [ExonData("1", NO_STRAND, 1),
                          ExonData("2", NO_STRAND, 2)]
        weights = calculate_node_weight(exon_data_list)
        self.assertTrue(weights == (3,0))
        # test combination
        exon_data_list = [ExonData("1", POS_STRAND, 10),
                          ExonData("2", NEG_STRAND, 20),
                          ExonData("3", NO_STRAND, 30)]
        weights = calculate_node_weight(exon_data_list)
        self.assertTrue(weights == (20,40))
        

#    def test_strand_degree(self):
#        """ensure that in_degree/out_degree are correctly ignoring nodes
#        on the opposite strand"""        
#        ig = TranscriptGraph()
#        t1 = make_transcript(((0, 100), (500,600), (900, 1000)), strand="+")
#        t2 = make_transcript(((300, 400), (500,600), (700,800)), strand="-")
#        ig.add_transcripts([t1, t2])
#        G = ig.G        
#        self.assertTrue(get_in_degree(G, Exon(0, 100), POS_STRAND) == 0)
#        self.assertTrue(get_out_degree(G, Exon(0, 100), POS_STRAND) == 1)
#        self.assertTrue(get_in_degree(G, Exon(500,600), POS_STRAND) == 1)
#        self.assertTrue(get_out_degree(G, Exon(500,600), POS_STRAND) == 1)
#        self.assertTrue(get_in_degree(G, Exon(900,1000), POS_STRAND) == 1)
#        self.assertTrue(get_out_degree(G, Exon(900,1000), POS_STRAND) == 0)
#        self.assertTrue(get_in_degree(G, Exon(300,400), NEG_STRAND) == 1)
#        self.assertTrue(get_out_degree(G, Exon(300,400), NEG_STRAND) == 0)
#        self.assertTrue(get_in_degree(G, Exon(500,600), NEG_STRAND) == 1)
#        self.assertTrue(get_out_degree(G, Exon(500,600), NEG_STRAND) == 1)
#        self.assertTrue(get_in_degree(G, Exon(700,800), NEG_STRAND) == 0)
#        self.assertTrue(get_out_degree(G, Exon(700,800), NEG_STRAND) == 1)

#    def test_assembly(self):
#        # exons do not overlap
#        test_basename = "assembly1"
#        gtf_file = test_basename + ".gtf"
#        dot_file = test_basename + ".dot"
#        isoform_graph = read_gtf(gtf_file)
#        write_dot(isoform_graph, "a")
#        for strand, gene_id, tss_id, score, path in isoform_graph.assemble(100):
#            print strand, gene_id, tss_id, score, path


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()