'''
Created on Feb 16, 2011

@author: mkiyer
'''
import logging
import unittest
import os
import networkx as nx
import matplotlib.pyplot as plt

from assemblyline.lib.base import Exon, POS_STRAND, NEG_STRAND, NO_STRAND
from assemblyline.lib.transcript_graph import TranscriptGraph
from assemblyline.lib.assembler import NODE_DENSITY, EDGE_OUT_FRAC, calculate_edge_attrs
from assemblyline.assemblyline import collapse_contiguous_nodes

from test_base import make_transcript, read_gtf, write_dot, get_dot_path

logging.basicConfig(level=logging.DEBUG,
                    #filename="/exds/users/mkiyer/projects/assemblyline/test/test.log",
                    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")

#class TestCollapseContiguousNodes(unittest.TestCase):
#    
#    def test_collapse1(self):
#        """
#        test function to merge genomically contiguous nodes
#        """
#        # plus strand
#        path = [Exon(0,10), Exon(10,15), Exon(16,20), Exon(20,30)]
#        newpath = collapse_contiguous_nodes(path, POS_STRAND)
#        correctpath = [Exon(0,15), Exon(16,30)]
#        for x,y in zip(newpath, correctpath):
#            self.assertTrue(x == y)
#        # minus strand
#        path = [Exon(50,60), Exon(49,50), Exon(48,49), Exon(40,45), 
#                Exon(39,40), Exon(0,10)]
#        newpath = collapse_contiguous_nodes(path, NEG_STRAND)
#        correctpath = [Exon(48,60), Exon(39,45), Exon(0,10)]
#        for x,y in zip(newpath, correctpath):
#            self.assertTrue(x == y)
    
#class TestCalculateWeights(unittest.TestCase):
#
#    def test_build_pos_strand_graph(self):
#        '''
#        test a simple graph with only forward strand transcripts
#        '''
#        tg = TranscriptGraph()
#        t1 = make_transcript(((0, 100), (200,300)), id="1", strand="+", score=20)
#        t2 = make_transcript(((0, 100), (500,600)), id="2", strand="+", score=30)
#        t3 = make_transcript(((0, 100), (700,800)), id="3", strand="+", score=50)
#        tg.add_transcripts([t1, t2, t3])
#        GG = build_strand_specific_graphs(tg.G)        
#        #for n,d in GG[POS_STRAND].nodes_iter(data=True):
#        #    print 'NODE', n,d
#        #for u,v,d in GG[POS_STRAND].edges_iter(data=True):
#        #    print u,v,d
#        calculate_edge_attrs(GG[POS_STRAND])
#        calculate_edge_attrs(GG[NEG_STRAND])        
#        # ensure no neg strand transcripts
#        self.assertTrue(len(GG[NEG_STRAND]) == 0)
#        G = GG[POS_STRAND]
#        # check node weights
#        self.assertAlmostEqual(G.node[t1.exons[0]][NODE_WEIGHT], 50)
#        self.assertAlmostEqual(G.node[t1.exons[1]][NODE_WEIGHT], 10)
#        self.assertAlmostEqual(G.node[t2.exons[1]][NODE_WEIGHT], 15)
#        self.assertAlmostEqual(G.node[t3.exons[1]][NODE_WEIGHT], 25)
#        # check edge weight computation
#        self.assertAlmostEqual(G.edge[t1.exons[0]][t1.exons[1]][EDGE_OUT_FRAC], 0.20)
#        self.assertAlmostEqual(G.edge[t1.exons[0]][t2.exons[1]][EDGE_OUT_FRAC], 0.30)
#        self.assertAlmostEqual(G.edge[t1.exons[0]][t3.exons[1]][EDGE_OUT_FRAC], 0.50)
#
#    def test_build_both_strand_graphs(self):
#        '''
#        test that transcript graph with both pos and neg strand transcripts
#        is property divided into forward strand and reverse strand graphs, 
#        and that node/edge weights are computed correctly
#        '''
#        tg = TranscriptGraph()
#        t1 = make_transcript(((0, 100), (200,300)), id="1", strand="+", score=20)
#        t2 = make_transcript(((0, 100), (500,600)), id="2", strand="+", score=30)
#        t3 = make_transcript(((0, 100), (700,800)), id="3", strand="+", score=50)
#        t4 = make_transcript(((0, 100), (200,300)), id="4", strand="-", score=20)
#        t5 = make_transcript(((0, 100), (500,600)), id="5", strand="-", score=30)
#        t6 = make_transcript(((0, 100), (700,800)), id="6", strand="-", score=50)
#        tg.add_transcripts([t1, t2, t3, t4, t5, t6])
#        GG = build_strand_specific_graphs(tg.G)
#        calculate_edge_attrs(GG[POS_STRAND])
#        calculate_edge_attrs(GG[NEG_STRAND])
#        # check number of nodes
#        self.assertTrue(len(GG[POS_STRAND]) == 4)
#        self.assertTrue(len(GG[NEG_STRAND]) == 4)
#        # check plus strand graph
#        G = GG[POS_STRAND]
#        self.assertAlmostEqual(G.node[t1.exons[0]][NODE_WEIGHT], 50)
#        self.assertAlmostEqual(G.node[t1.exons[1]][NODE_WEIGHT], 10)
#        self.assertAlmostEqual(G.node[t2.exons[1]][NODE_WEIGHT], 15)
#        self.assertAlmostEqual(G.node[t3.exons[1]][NODE_WEIGHT], 25)
#        self.assertAlmostEqual(G.edge[t1.exons[0]][t1.exons[1]][EDGE_OUT_FRAC], 0.20)
#        self.assertAlmostEqual(G.edge[t1.exons[0]][t2.exons[1]][EDGE_OUT_FRAC], 0.30)
#        self.assertAlmostEqual(G.edge[t1.exons[0]][t3.exons[1]][EDGE_OUT_FRAC], 0.50)
#        # check minus strand graph
#        G = GG[NEG_STRAND]
#        self.assertAlmostEqual(G.node[t1.exons[0]][NODE_WEIGHT], 50)
#        self.assertAlmostEqual(G.node[t1.exons[1]][NODE_WEIGHT], 10)
#        self.assertAlmostEqual(G.node[t2.exons[1]][NODE_WEIGHT], 15)
#        self.assertAlmostEqual(G.node[t3.exons[1]][NODE_WEIGHT], 25)
#        self.assertAlmostEqual(G.edge[t1.exons[1]][t1.exons[0]][EDGE_OUT_FRAC], 1.0)
#        self.assertAlmostEqual(G.edge[t2.exons[1]][t1.exons[0]][EDGE_OUT_FRAC], 1.0)
#        self.assertAlmostEqual(G.edge[t3.exons[1]][t1.exons[0]][EDGE_OUT_FRAC], 1.0)
#
#    def test_build_with_unstranded_transcripts(self):
#        '''
#        test that transcript graph with some unstranded transcripts
#        is property divided into forward strand and reverse strand graphs, 
#        and that node/edge weights are computed correctly
#        '''
#        tg = TranscriptGraph()
#        t1 = make_transcript(((0, 100), (200,300)), id="1", strand="+", score=300)
#        t2 = make_transcript(((600, 700), (900,1000)), id="2", strand="-", score=100)
#        t3 = make_transcript(((0, 1000),), id="3", strand=".", score=400)
#        tg.add_transcripts([t1, t2, t3])        
#        GG = build_strand_specific_graphs(tg.G)
#        calculate_edge_attrs(GG[POS_STRAND])
#        calculate_edge_attrs(GG[NEG_STRAND])
#        # check number of nodes
#        self.assertTrue(len(GG[POS_STRAND]) == 7)
#        self.assertTrue(len(GG[NEG_STRAND]) == 7)
#        # the unstranded transcript is split 75/25 between pos/neg strands due to
#        # the coverage in both         
#        pos_frac = 300.0/(300.0+100.0)
#        # check graph
#        GP = GG[POS_STRAND]
#        GN = GG[NEG_STRAND]
#        # pos exons
#        self.assertAlmostEqual(GP.node[t1.exons[0]][NODE_WEIGHT], (pos_frac*0.1*400) + (0.5*300))
#        self.assertAlmostEqual(GP.node[t1.exons[1]][NODE_WEIGHT], (pos_frac*0.1*400) + (0.5*300))    
#        # neg exons
#        self.assertAlmostEqual(GN.node[t2.exons[0]][NODE_WEIGHT], ((1-pos_frac)*0.1*400) + (0.5*100))
#        self.assertAlmostEqual(GN.node[t2.exons[1]][NODE_WEIGHT], ((1-pos_frac)*0.1*400) + (0.5*100))
#        # unstranded places    
#        self.assertAlmostEqual(GP.node[Exon(100,200)][NODE_WEIGHT], pos_frac*0.1*400)
#        self.assertAlmostEqual(GN.node[Exon(100,200)][NODE_WEIGHT], (1-pos_frac)*0.1*400)
#        # edge weights
#        self.assertAlmostEqual(GP.edge[t1.exons[0]][t1.exons[1]][EDGE_OUT_FRAC], 5./6.)
#        self.assertAlmostEqual(GP.edge[t1.exons[0]][Exon(100,200)][EDGE_OUT_FRAC], 1./6.)
#        self.assertAlmostEqual(GN.edge[t2.exons[1]][t2.exons[0]][EDGE_OUT_FRAC], 5./6.)
#        self.assertAlmostEqual(GN.edge[t2.exons[1]][Exon(700,900)][EDGE_OUT_FRAC], 1./6.)


#class TestAssembler(unittest.TestCase):
#    
#    def test_find_paths(self):
#        '''
#        test that assembler finds the intended "best" paths
#        '''
#        tg = TranscriptGraph()
#        t1 = make_transcript(((0, 100), (200,300), (900,1000)), id="1", strand="+", score=900)
#        t2 = make_transcript(((0, 100), (500,600), (900,1000)), id="2", strand="+", score=600)
#        t3 = make_transcript(((0, 100), (700,800), (900,1000)), id="3", strand="+", score=300)
#        tg.add_transcripts([t1, t2, t3])
#        results = list(tg.assemble(max_paths=10, fraction_major_path=0))
#        self.assertTrue(len(results) == 3)
#        #print len(results)
#        #yield strand, gene_id, tss_id, score, path
#
#    def test_actual_data1(self):
#        """test an example taken from real data"""
#        # exons do not overlap
#        test_basename = "assembly2"
#        gtf_file = test_basename + ".gtf"
#        dot_file = test_basename + ".dot"
#        txgraph = read_gtf(gtf_file)
#        write_dot(txgraph, "a")
#        from ..assembler import build_strand_specific_graphs
#        GG = build_strand_specific_graphs(txgraph.G)
#        nx.write_dot(GG[0], get_dot_path("b"))
#        nx.write_dot(GG[1], get_dot_path("c"))


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()