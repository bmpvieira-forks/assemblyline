'''
Created on Nov 3, 2010

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


class TestMergeTranscripts(unittest.TestCase):

    def testMergeIdentical(self):
        gtf_file = os.path.join(os.path.dirname(__file__), "merge_identical.gtf")
        for locus_transcripts in parse_gtf(open(gtf_file)):
            isoform_graph = IsoformGraph.from_transcripts(locus_transcripts)
            isoform_graph.collapse()
            Gcorrect = nx.read_dot(os.path.join(os.path.dirname(__file__), "merge_identical.dot"))            
            self.assertTrue(all(str(n) in Gcorrect for n in isoform_graph.G))
            self.assertTrue(nx.is_isomorphic(isoform_graph.G, Gcorrect)) 

    def testMergeSeparateStrands(self):
        gtf_file = os.path.join(os.path.dirname(__file__), "merge_separate_strands.gtf")
        for locus_transcripts in parse_gtf(open(gtf_file)):
            isoform_graph = IsoformGraph.from_transcripts(locus_transcripts) 
            isoform_graph.collapse()
            Gcorrect = nx.read_dot(os.path.join(os.path.dirname(__file__), "merge_separate_strands.dot"))            
            self.assertTrue(all(str(n) in Gcorrect for n in isoform_graph.G))
            self.assertTrue(nx.is_isomorphic(isoform_graph.G, Gcorrect)) 

    def testMergeSeparateChromosomes(self):
        gtf_file = os.path.join(os.path.dirname(__file__), "merge_separate_chroms.gtf")
        for locus_transcripts in parse_gtf(open(gtf_file)):
            isoform_graph = IsoformGraph.from_transcripts(locus_transcripts) 
            isoform_graph.collapse()
            self.assertTrue(isoform_graph.G.number_of_nodes() == 1)
    
    def testMergePath(self):
        gtf_file = os.path.join(os.path.dirname(__file__), "merge_path.gtf")
        for locus_transcripts in parse_gtf(open(gtf_file)):
            isoform_graph = IsoformGraph.from_transcripts(locus_transcripts) 
            isoform_graph.collapse()
            Gcorrect = nx.read_dot(os.path.join(os.path.dirname(__file__), "merge_path.dot"))
            self.assertTrue(all(str(n) in Gcorrect for n in isoform_graph.G))
            self.assertTrue(nx.is_isomorphic(isoform_graph.G, Gcorrect))

    def testMergePath2(self):
        gtf_file = os.path.join(os.path.dirname(__file__), "merge_path2.gtf")
        for locus_transcripts in parse_gtf(open(gtf_file)):
            isoform_graph = IsoformGraph.from_transcripts(locus_transcripts) 
            isoform_graph.collapse()
            Gcorrect = nx.read_dot(os.path.join(os.path.dirname(__file__), "merge_path2.dot"))
            self.assertTrue(all(str(n) in Gcorrect for n in isoform_graph.G))
            self.assertTrue(nx.is_isomorphic(isoform_graph.G, Gcorrect))            

    def testMergePath3(self):
        """test merging exons that are touching but not overlapping"""
        gtf_file = os.path.join(os.path.dirname(__file__), "merge_path3.gtf")
        for locus_transcripts in parse_gtf(open(gtf_file)):
            isoform_graph = IsoformGraph.from_transcripts(locus_transcripts) 
            isoform_graph.collapse()
            Gcorrect = nx.read_dot(os.path.join(os.path.dirname(__file__), "merge_path3.dot"))
            self.assertTrue(all(str(n) in Gcorrect for n in isoform_graph.G))
            self.assertTrue(nx.is_isomorphic(isoform_graph.G, Gcorrect))

    def testMergePath4(self):
        """test merging two transcripts, where one is a truncated version of the other"""
        gtf_file = os.path.join(os.path.dirname(__file__), "merge_path4.gtf")
        for locus_transcripts in parse_gtf(open(gtf_file)):
            isoform_graph = IsoformGraph.from_transcripts(locus_transcripts)
            isoform_graph.collapse()
            Gcorrect = nx.read_dot(os.path.join(os.path.dirname(__file__), "merge_path4.dot"))
            self.assertTrue(all(str(n) in Gcorrect for n in isoform_graph.G))
            self.assertTrue(nx.is_isomorphic(isoform_graph.G, Gcorrect))

    def testMergePath5(self):
        """test merging two transcripts, where one is a truncated version of the other"""
        gtf_file = os.path.join(os.path.dirname(__file__), "merge_path5.gtf")
        for locus_transcripts in parse_gtf(open(gtf_file)):
            isoform_graph = IsoformGraph.from_transcripts(locus_transcripts)
            isoform_graph.collapse()
            Gcorrect = nx.read_dot(os.path.join(os.path.dirname(__file__), "merge_path5.dot"))
            self.assertTrue(all(str(n) in Gcorrect for n in isoform_graph.G))
            self.assertTrue(nx.is_isomorphic(isoform_graph.G, Gcorrect))

    def testMergePath6(self):
        """test a complex locus on chromosome 1 that gave errors due to problems with threading"""
        gtf_file = os.path.join(os.path.dirname(__file__), "merge_path6.gtf")
        for locus_transcripts in parse_gtf(open(gtf_file)):
            isoform_graph = IsoformGraph.from_transcripts(locus_transcripts)
            isoform_graph.collapse(trim=True, overhang_threshold=15)
            #nx.write_dot(isoform_graph.G, "merge_path6.dot")
            Gcorrect = nx.read_dot(os.path.join(os.path.dirname(__file__), "merge_path6.dot"))
            self.assertTrue(all(str(n) in Gcorrect for n in isoform_graph.G))
            self.assertTrue(nx.is_isomorphic(isoform_graph.G, Gcorrect))
            # TODO: need to test isoform generation
            break

    def testMergePath7(self):
        """test a complex locus on chromosome 1 that gave errors due to problems with threading"""
        gtf_file = os.path.join(os.path.dirname(__file__), "merge_path7.gtf")
        for locus_transcripts in parse_gtf(open(gtf_file)):
            isoform_graph = IsoformGraph.from_transcripts(locus_transcripts)
            nx.spring_layout(isoform_graph.G)
            nx.draw(isoform_graph.G)
            plt.show()                   
            isoform_graph.collapse(trim=True, overhang_threshold=15)

            nx.spring_layout(isoform_graph.G)
            nx.draw(isoform_graph.G)
            plt.show()
            nx.write_dot(isoform_graph.G, "tmp.dot")
            Gcorrect = nx.read_dot(os.path.join(os.path.dirname(__file__), "merge_path7.dot"))
            self.assertTrue(all(str(n) in Gcorrect for n in isoform_graph.G))
            self.assertTrue(nx.is_isomorphic(isoform_graph.G, Gcorrect))
            # TODO: need to test isoform generation
#
    def testMergeWithStrandlessExons(self):
        """test merging two stranded transcripts with an unstranded transcript"""
        gtf_file = os.path.join(os.path.dirname(__file__), "merge_nostrand.gtf")
        for locus_transcripts in parse_gtf(open(gtf_file)):
            isoform_graph = IsoformGraph.from_transcripts(locus_transcripts) 
            isoform_graph.collapse()
            Gcorrect = nx.read_dot(os.path.join(os.path.dirname(__file__), "merge_nostrand.dot"))
            self.assertTrue(all(str(n) in Gcorrect for n in isoform_graph.G))
            self.assertTrue(nx.is_isomorphic(isoform_graph.G, Gcorrect))

    def testMultipleMerges(self):
        """test merge involving multiple paths"""
        gtf_file = os.path.join(os.path.dirname(__file__), "merge_multi.gtf")
        for locus_transcripts in parse_gtf(open(gtf_file)):
            isoform_graph = IsoformGraph.from_transcripts(locus_transcripts) 
            isoform_graph.collapse()            
            Gcorrect = nx.read_dot(os.path.join(os.path.dirname(__file__), "merge_multi.dot"))
            self.assertTrue(all(str(n) in Gcorrect for n in isoform_graph.G))
            self.assertTrue(nx.is_isomorphic(isoform_graph.G, Gcorrect))

    def testMultipleMerges2(self):
        """test merge involving multiple paths"""
        gtf_file = os.path.join(os.path.dirname(__file__), "merge_multi2.gtf")
        for locus_transcripts in parse_gtf(open(gtf_file)):
            isoform_graph = IsoformGraph.from_transcripts(locus_transcripts) 
            isoform_graph.collapse()
            Gcorrect = nx.read_dot(os.path.join(os.path.dirname(__file__), "merge_multi2.dot"))
            self.assertTrue(all(str(n) in Gcorrect for n in isoform_graph.G))
            self.assertTrue(nx.is_isomorphic(isoform_graph.G, Gcorrect))

    def testIncompatibleMerge(self):
        """test merge involving intron-incompatible paths"""
        gtf_file = os.path.join(os.path.dirname(__file__), "merge_incompatible.gtf")
        for locus_transcripts in parse_gtf(open(gtf_file)):
            isoform_graph = IsoformGraph.from_transcripts(locus_transcripts) 
            isoform_graph.collapse()
            nx.write_dot(isoform_graph.G, "merge_incompatible.dot")            
            Gcorrect = nx.read_dot(os.path.join(os.path.dirname(__file__), "merge_incompatible.dot"))
            self.assertTrue(all(str(n) in Gcorrect for n in isoform_graph.G))
            self.assertTrue(nx.is_isomorphic(isoform_graph.G, Gcorrect))

    def testMergeComplexPath(self):
        """test merging transcripts where one transcript encompasses the other"""
        gtf_file = os.path.join(os.path.dirname(__file__), "merge_complex_path.gtf")
        for locus_transcripts in parse_gtf(open(gtf_file)):
            isoform_graph = IsoformGraph.from_transcripts(locus_transcripts)
            isoform_graph.collapse()
            Gcorrect = nx.read_dot(os.path.join(os.path.dirname(__file__), "merge_complex_path.dot"))
            self.assertTrue(all(str(n) in Gcorrect for n in isoform_graph.G))
            self.assertTrue(nx.is_isomorphic(isoform_graph.G, Gcorrect))

    def testTrim(self):
        """test trimming intron-incompatible paths"""
        gtf_file = os.path.join(os.path.dirname(__file__), "trim.gtf")
        for locus_transcripts in parse_gtf(open(gtf_file)):
            isoform_graph = IsoformGraph.from_transcripts(locus_transcripts)

            nx.spring_layout(isoform_graph.G)
            nx.draw(isoform_graph.G)
            plt.show()

            isoform_graph.collapse(trim=True, overhang_threshold=15)

            nx.spring_layout(isoform_graph.G)
            nx.draw(isoform_graph.G)
            plt.show()
            
            #Gcorrect = nx.read_dot(os.path.join(os.path.dirname(__file__), "trim.dot"))
            #self.assertTrue(all(str(n) in Gcorrect for n in isoform_graph.G))
            #self.assertTrue(nx.is_isomorphic(isoform_graph.G, Gcorrect))

#    def testThreading(self):
#        """test threading transcripts to find paths through graph"""
#        gtf_file = os.path.join(os.path.dirname(__file__), "thread_transcripts.gtf")
#        for locus_transcripts in parse_gtf(open(gtf_file)):
#            isoform_graph = IsoformGraph.from_transcripts(locus_transcripts)
#            isoform_graph.get_isoforms(locus_transcripts)
            # TODO: do a test!
            #nx.write_dot(isoform_graph.G, "merge_path7.dot")

            #nx.spectral_layout(isoform_graph.G)
            #nx.spring_layout(isoform_graph.G)
            #nx.draw(isoform_graph.G)
            #plt.show()
#            nx.write_dot(isoform_graph.G, "merge_path6.dot")

#            nx.write_dot(isoform_graph.G, "merge_path5.dot")            
#            nx.write_dot(isoform_graph.G, "merge_complex_path.dot")            
#            nx.spring_layout(isoform_graph.G)
#            nx.draw(isoform_graph.G)
#            plt.show()
#            nx.write_dot(isoform_graph.G, "merge_incompatible.dot")            
#            nx.spring_layout(isoform_graph.G)
#            nx.draw(isoform_graph.G)
#            plt.show()
            #nx.write_dot(isoform_graph.G, "merge_multi2.dot")            
            #Gcorrect = nx.read_dot(os.path.join(os.path.dirname(__file__), "merge_nostrand.dot"))
            #nx.write_dot(isoform_graph.G, "merge_path3.dot")
            #self.assertTrue(nx.is_isomorphic(isoform_graph.G, Gcorrect)) 
#            nx.draw(isoform_graph.G)
#            plt.show()
            #nx.spectral_layout(isoform_graph.G)
            #nx.spring_layout(isoform_graph.G)
            #nx.draw(isoform_graph.G)
            #plt.show()
            #Gcorrect = nx.read_dot(os.path.join(os.path.dirname(__file__), "trim.dot"))
            #self.assertTrue(all(str(n) in Gcorrect for n in isoform_graph.G))
            #self.assertTrue(nx.is_isomorphic(isoform_graph.G, Gcorrect))
#            nx.spring_layout(isoform_graph.G)
#            nx.draw(isoform_graph.G)
#            plt.show()
#            nx.write_dot(isoform_graph.G, "merge_incompatible.dot")            
#            nx.spring_layout(isoform_graph.G)
#            nx.draw(isoform_graph.G)
#            plt.show()
            #nx.write_dot(isoform_graph.G, "merge_multi2.dot")            
            #Gcorrect = nx.read_dot(os.path.join(os.path.dirname(__file__), "merge_nostrand.dot"))
            #nx.write_dot(isoform_graph.G, "merge_path3.dot")
            #self.assertTrue(nx.is_isomorphic(isoform_graph.G, Gcorrect)) 
#            nx.draw(isoform_graph.G)
#            plt.show()


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()