'''
Created on Mar 13, 2011

@author: mkiyer
'''
import logging
import unittest
import os
import networkx as nx
import matplotlib.pyplot as plt

logging.basicConfig(level=logging.DEBUG,
                    #filename="/exds/users/mkiyer/projects/assemblyline/test/test.log",
                    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")

from assemblyline.lib.transcript import Exon, POS_STRAND, NEG_STRAND, NO_STRAND
from assemblyline.lib.assembler import calculate_edge_attrs, smoothen_graph_leastsq

def example1():
    G = nx.DiGraph()
    nodes = [Exon(0, 100),
             Exon(200, 300),
             Exon(400, 500),
             Exon(525, 575),
             Exon(600, 700),
             Exon(800, 1000),
             Exon(1100,1200)]
    G.add_node(nodes[0], weight=100)
    G.add_node(nodes[1], weight=100)
    G.add_node(nodes[2], weight=200)
    G.add_node(nodes[3], weight=400)
    G.add_node(nodes[4], weight=300)
    G.add_node(nodes[5], weight=10)
    G.add_node(nodes[6], weight=10)
    G.add_edge(nodes[0], nodes[2])
    G.add_edge(nodes[1], nodes[2])
    G.add_edge(nodes[2], nodes[3])
    G.add_edge(nodes[3], nodes[4])
    G.add_edge(nodes[4], nodes[5])
    G.add_edge(nodes[4], nodes[6])
    return G

def example2():
    G = nx.DiGraph()
    # 1
    G.add_node(Exon(0, 100), weight=100)
    # 2
    G.add_node(Exon(200, 300), weight=90)
    G.add_edge(Exon(0, 100), Exon(200, 300), density=75)
    # 3
    G.add_node(Exon(400, 500), weight=50)
    G.add_edge(Exon(0, 100), Exon(400, 500), density=25)
    # 4
    G.add_node(Exon(600, 700), weight=200)
    G.add_edge(Exon(200, 300), Exon(600, 700), density=100)
    G.add_edge(Exon(400, 500), Exon(600, 700), density=50)
    # 5
    G.add_node(Exon(800, 900), weight=150)
    G.add_edge(Exon(400, 500), Exon(800, 900), density=50)
    # 6
    G.add_node(Exon(1000, 1100), weight=500)
    G.add_edge(Exon(600, 700), Exon(1000, 1100), density=100)
    G.add_edge(Exon(800, 900), Exon(1000, 1100), density=100)
    # 7
    G.add_node(Exon(-200, -100), weight=100)
    G.add_edge(Exon(-200, -100), Exon(400,500), density=100)
    return G

class TestSmoothen(unittest.TestCase):

    def test_smoothen_leastsq(self):
        G = example2()
        calculate_edge_attrs(G)
        smoothen_graph_leastsq(G)
        # check result
        self.assertEqual(int(G.node[Exon(0, 100)]['weight']), 182)
        self.assertEqual(int(G.node[Exon(200, 300)]['weight']), 136)
        self.assertEqual(int(G.node[Exon(400, 500)]['weight']), 190)
        self.assertEqual(int(G.node[Exon(600, 700)]['weight']), 231)
        self.assertEqual(int(G.node[Exon(800, 900)]['weight']), 95)
        self.assertEqual(int(G.node[Exon(1000, 1100)]['weight']), 326)
        self.assertEqual(int(G.node[Exon(-200, -100)]['weight']), 144)
        return G
        
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()