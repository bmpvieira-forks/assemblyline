'''
Created on Dec 13, 2011

@author: mkiyer
'''
import unittest

import networkx as nx
import numpy as np

from assemblyline.lib.transcript import Exon
from assemblyline.lib.assembler_base import NODE_DENSITY, NODE_LENGTH, EDGE_OUT_FRAC, EDGE_IN_FRAC, \
    MIN_NODE_DENSITY
from assemblyline.lib.path_finder import calculate_edge_fractions, \
    dynprog_search, traceback, find_path, init_tmp_attrs, get_seed_subgraph, PATH_DICT
    
def example1():
    G = nx.DiGraph()
    # add nodes
    G.add_node('SOURCE', attr_dict={NODE_DENSITY: MIN_NODE_DENSITY,
                                    NODE_LENGTH: 0})
    G.add_node('A', attr_dict={NODE_DENSITY: 20, NODE_LENGTH: 1})
    G.add_node('B', attr_dict={NODE_DENSITY: 10, NODE_LENGTH: 1})
    G.add_node('C', attr_dict={NODE_DENSITY: 10, NODE_LENGTH: 1})
    G.add_node('D', attr_dict={NODE_DENSITY: 30, NODE_LENGTH: 1})
    G.add_node('H', attr_dict={NODE_DENSITY: 20, NODE_LENGTH: 1})
    G.add_node('E', attr_dict={NODE_DENSITY: 10, NODE_LENGTH: 1})
    G.add_node('F', attr_dict={NODE_DENSITY: 20, NODE_LENGTH: 1})
    G.add_node('G', attr_dict={NODE_DENSITY: 40, NODE_LENGTH: 1})
    G.add_node('SINK', attr_dict={NODE_DENSITY: MIN_NODE_DENSITY, 
                                  NODE_LENGTH: 0})
    # add edges
    G.add_edges_from([('SOURCE','A'), ('A','B'), ('A','C'), ('B','H'), 
                      ('B','D'), ('C','D'), ('H','E'), ('D','E'), 
                      ('D','F'), ('E','G'), ('F','G'), ('G','SINK')])
    return G

def example2():
    G = nx.DiGraph()
    # add nodes
    G.add_node('SOURCE', attr_dict={NODE_DENSITY: MIN_NODE_DENSITY,
                                    NODE_LENGTH: 0})
    G.add_node('A', attr_dict={NODE_DENSITY: 30, NODE_LENGTH: 1})
    G.add_node('B', attr_dict={NODE_DENSITY: 20, NODE_LENGTH: 1})
    G.add_node('C', attr_dict={NODE_DENSITY: 15, NODE_LENGTH: 1})
    G.add_node('D', attr_dict={NODE_DENSITY: 5, NODE_LENGTH: 1})
    G.add_node('E', attr_dict={NODE_DENSITY: 20, NODE_LENGTH: 1})
    G.add_node('F', attr_dict={NODE_DENSITY: 20, NODE_LENGTH: 1})
    G.add_node('G', attr_dict={NODE_DENSITY: 10, NODE_LENGTH: 1})
    G.add_node('SINK', attr_dict={NODE_DENSITY: MIN_NODE_DENSITY, 
                                  NODE_LENGTH: 0})
    # add edges
    G.add_edges_from([('SOURCE','A'), ('A','B'), ('A','G'), ('B','C'), 
                      ('B','D'), ('C', 'E'), ('D', 'E'), ('E', 'F'),
                      ('F', 'SINK'), ('G', 'F')]) 
    return G
    
    
class TestPathFinder(unittest.TestCase):

    def test_calculate_edge_fractions(self):
        G = example1()
        calculate_edge_fractions(G)
        # test some random edges
        self.assertAlmostEqual(G['B']['D'][EDGE_OUT_FRAC], 3./5.)
        self.assertAlmostEqual(G['D']['F'][EDGE_OUT_FRAC], 2./3.)
        self.assertAlmostEqual(G['C']['D'][EDGE_IN_FRAC], 1./2.)
        self.assertAlmostEqual(G['E']['G'][EDGE_OUT_FRAC], 1./1.)
        self.assertAlmostEqual(G['E']['G'][EDGE_IN_FRAC], 1./3.)
        self.assertAlmostEqual(G['F']['G'][EDGE_IN_FRAC], 2./3.)

    def test_find_path(self):
        G = example1()
        calculate_edge_fractions(G)
        init_tmp_attrs(G)
        dynprog_search(G, 'SOURCE')
        # hard-code the correct result
        PATH_DICTS = {'A': {'SOURCE': (0.0, 0)},
                      'C': {'A': (10.0, 1.0)},
                      'B': {'A': (10.0, 1.0)},
                      'E': {'H': (34.0, 3.0), 'D': (25.0, 3.0)},
                      'D': {'C': (20.0, 2.0), 'B': (16.0, 2.0)},
                      'G': {'E': (38.0, 4.0), 'F': (50.0, 4.0)},
                      'F': {'D': (30.0, 3.0)},
                      'H': {'B': (14.0, 2.0)},
                      'SOURCE': {},
                      'SINK': {'G': (76.66666666666666, 5.0)}}
        for n,d in G.nodes_iter(data=True):
            self.assertAlmostEqual(d[PATH_DICT], PATH_DICTS[n])
        path, weight, length = traceback(G, 'SINK', inclusive=True)
        self.assertEqual(path, ('SOURCE', 'A', 'C', 'D', 'F', 'G', 'SINK'))
        self.assertAlmostEqual(weight, (76 + 2./3.))
        self.assertEqual(length, 5)
    
    def test_find_path_reverse(self):
        G = example1()
        calculate_edge_fractions(G)
        init_tmp_attrs(G)
        # try to find the reverse path
        path, weight, length = find_path(G, 'SOURCE', 'SINK', reverse=True)
        self.assertEqual(path, ('SINK', 'G', 'F', 'D', 'C', 'A', 'SOURCE'))
        self.assertAlmostEqual(weight, (76 + 2./3.))
        self.assertEqual(length, 5)

    def test_get_seed_subgraph(self):
        G = example2()
        calculate_edge_fractions(G)
        init_tmp_attrs(G)
        A = (np.identity(len(G)) - nx.adj_matrix(G)).I
        nodes = G.nodes()
        node_indexes = dict((n,i) for i,n in enumerate(nodes))        
        H = get_seed_subgraph(G, A, nodes, node_indexes, 'C')
        # check subgraph
        for n in ('SOURCE', 'A', 'B', 'C', 'E', 'F', 'SINK'):            
            self.assertTrue(n in H)
        for n in ('D', 'G'):
            self.assertFalse(n in H)
        
        



if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()