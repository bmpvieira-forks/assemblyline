'''
Created on Dec 7, 2011

@author: mkiyer
'''
import unittest

import networkx as nx

from assemblyline.lib.transcript import Exon
from assemblyline.lib.collapse_chains import can_collapse_contiguous

class TestCollapse(unittest.TestCase):
        
    def test_collapse_contiguous(self):
        # make graph
        G = nx.DiGraph()
        a = Exon(0,500)
        b = Exon(500,600)
        c = Exon(600,700)
        d = Exon(700,800)
        e = Exon(800,900)
        f = Exon(1500,1600)
        G.add_path([a,b,c,d,e,f])
        self.assertTrue(can_collapse_contiguous(G,a,b))
        self.assertTrue(can_collapse_contiguous(G,b,c))
        self.assertTrue(can_collapse_contiguous(G,c,d))
        self.assertTrue(can_collapse_contiguous(G,d,e))
        self.assertFalse(can_collapse_contiguous(G,e,f))
        # add edge to invalidate collapsing
        G.add_edge(a,d)
        self.assertFalse(can_collapse_contiguous(G,a,b))
        self.assertTrue(can_collapse_contiguous(G,b,c))
        self.assertFalse(can_collapse_contiguous(G,c,d))
        self.assertTrue(can_collapse_contiguous(G,d,e))
        # add reverse edge which shouldn't affect anything
        G.add_edge(d,a)
        self.assertFalse(can_collapse_contiguous(G,a,b))
        self.assertTrue(can_collapse_contiguous(G,b,c))
        self.assertFalse(can_collapse_contiguous(G,c,d))
        self.assertTrue(can_collapse_contiguous(G,d,e))
        # add another edge to invalidate
        G.add_edge(f,b)
        self.assertFalse(can_collapse_contiguous(G,a,b))
        self.assertFalse(can_collapse_contiguous(G,b,c))
        self.assertFalse(can_collapse_contiguous(G,c,d))
        self.assertTrue(can_collapse_contiguous(G,d,e))

    


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()