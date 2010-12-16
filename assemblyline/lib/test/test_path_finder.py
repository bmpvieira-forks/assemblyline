'''
Created on Dec 8, 2010

@author: mkiyer
'''
import unittest
import os
import networkx as nx
import matplotlib.pyplot as plt

from isoform_graph import IsoformGraph
from transcript_parser import parse_gtf
from path_finder import get_isoforms
from cNode import Node

class Test(unittest.TestCase):

    def testFindPaths1(self):
        gtf_file = os.path.join(os.path.dirname(__file__), "findpaths1.gtf")
        for locus_transcripts in parse_gtf(open(gtf_file)):
            isoform_graph = IsoformGraph.from_transcripts(locus_transcripts)

            tx_id = 0
            for gene_id, tss_id, score, path in get_isoforms(isoform_graph.G, locus_transcripts,
                                                             fraction_major_path=0.000001):
                print gene_id, tss_id, tx_id, score, path
                tx_id += 1
            return
        
            id_score_map = get_transcript_score_map(locus_transcripts)
            H = transform_graph(isoform_graph.G, id_score_map)            
            start_nodes = [n for (n,d) in H.in_degree_iter()
                           if (d == 0)]       
            for start_node in start_nodes:
                for p in dfs_paths(H, start_node):
                    print 'final path', str(p)

            #nx.write_dot(isoform_graph.G, "merge_path7.dot")
            #nx.spectral_layout(isoform_graph.G)
            #nx.spring_layout(isoform_graph.G)
            #nx.draw(isoform_graph.G)
            #plt.show()



if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()