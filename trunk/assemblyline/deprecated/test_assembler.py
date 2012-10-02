'''
Created on Dec 18, 2011

@author: mkiyer
'''
import unittest
import networkx as nx

from assemblyline.lib.assemble.assembler import get_partial_path_kmers, add_dummy_start_end_nodes
from assemblyline.lib.transcript import Exon

class TestAssembler(unittest.TestCase):

    def test_enumerate_paths(self):
        G = nx.DiGraph()
        G.add_path([-1, 1, 2, 3, 4, 5, 6, -5, -6, -7])
        G.add_path([-1, 1, 3, 4, 7, 8, -2, -3, -4])
        paths = get_partial_path_kmers(G, (4,7), 4)
        for p in paths:
            print p

    def test_add_dummy_nodes(self):
        G = nx.DiGraph()
        G.add_path([1, 2, 3, 4, 5, 6])
        G.add_path([1, 3, 4, 7, 8])
        add_dummy_start_end_nodes(G, 3)
        all_edges = set(G.edges())
        correct_edges = set([(1, 2), (1, 3), (2, 3), (3, 4), (4, 5), 
                             (4, 7), (5, 6), (7, 8), 
                             (6, Exon(-3,-3)),
                             (8, Exon(-5,-5)),
                             (Exon(-2,-2), 1), 
                             (Exon(-5,-5), Exon(-6,-6)), 
                             (Exon(-3,-3), Exon(-4,-4)), 
                             (Exon(-1,-1), Exon(-2,-2))])
        self.assertEqual(all_edges, correct_edges)
        



#    def test_kmer_dicts(self):
#        partial_paths = ['ABCDEFG', 'ABCEFG', 'DEFG']
#        d = build_kmer_dicts(partial_paths, kmax=5)
#
#        for x in d:        
#            for n,kmer_dict in x.iteritems():
#                print n, sorted(kmer_dict.items())
        #for k,dicts in d.iteritems():
        #    print k, 'fwd', sorted(dicts[0].items())
        #    print k, 'rev', sorted(dicts[1].items())


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()