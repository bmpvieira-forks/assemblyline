'''
Created on Dec 6, 2010

@author: mkiyer
'''
import unittest
import os
import networkx as nx
import matplotlib.pyplot as plt
import time

from ..kclique import kcliques

import logging
logging.basicConfig(level=logging.DEBUG,
                    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")


class TestKClique(unittest.TestCase):

    def testBasic(self):
        G = nx.Graph()
        # make a 4-clique
        G.add_path([1,2,3,4,1])
        G.add_edges_from([(1,3),(2,4)])
        # find 2-cliques
        cliques = kcliques(G, 2)
        correct = set([frozenset([2, 4]), frozenset([1, 2]), frozenset([1, 4]), frozenset([3, 4]), frozenset([2, 3]), frozenset([1, 3])])
        self.assertTrue(cliques == correct)
        # find 3-cliques
        cliques = kcliques(G, 3)
        correct = set([frozenset([1, 3, 4]), frozenset([1, 2, 3]), frozenset([2, 3, 4]), frozenset([1, 2, 4])])
        self.assertTrue(cliques == correct)
        # find 4-cliques
        cliques = kcliques(G, 4)
        correct = set([frozenset([1, 2, 3, 4])]) 
        self.assertTrue(cliques == correct)
        # find 5-cliques (none exist)
        cliques = kcliques(G, 5)
        correct = set([frozenset([1, 2, 3, 4])]) 
        self.assertTrue(cliques == correct)        

    def testTriangles(self):
        G = nx.Graph()
        # create a large graph with many nodes and one triangle
        G.add_edges_from(list((i,i+1) for i in xrange(1,100)))
        G.add_edge(1,3)
        cliques = kcliques(G, 3)
        self.assertTrue(frozenset([1,2,3]) in cliques)
        self.assertTrue(len([c for c in cliques if len(c) > 2]) == 1)
        
    def testComplete(self):
        for n in xrange(2,10):
            G = nx.complete_graph(n)
            start = time.time()
            cliques = kcliques(G, n)
            end = time.time()
            logging.debug("k=%d time=%f" % (n, end-start))
            self.assertTrue(frozenset(range(n)) in cliques)
        
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()