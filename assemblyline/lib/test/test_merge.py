'''
Created on Nov 23, 2010

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
def get_path(filename):
    return os.path.join(os.path.dirname(__file__), DOT_DIR, filename)

def convert_ids_to_string(ids):
    return ','.join([str(id) for id in sorted(ids)])
    
def convert_attrs_to_strings(H):
    G = H.copy()
    for n in G:
        G.node[n]['ids'] = convert_ids_to_string(G.node[n]['ids'])
    return G

def write_dot(isoform_graph, filename):
    G = convert_attrs_to_strings(isoform_graph.G)
    nx.write_dot(G, get_path(filename))

class TestMergeTranscripts(unittest.TestCase):

    def compare_dot(self, isoform_graph, filename):
        Gcorrect = nx.read_dot(get_path(filename))
        G = convert_attrs_to_strings(isoform_graph.G)
        for n in G:
            node_str = str(n)
            self.assertTrue(node_str in Gcorrect)
            correct_ids = Gcorrect.node[node_str]['ids']
            test_ids = '"%s"' % (G.node[n]['ids'])
            self.assertTrue(correct_ids == test_ids)
        self.assertTrue(nx.is_isomorphic(isoform_graph.G, Gcorrect)) 

    def testIntronCompatibility(self):
        """ensure internal nodes are not marked as intron compatible"""
        gtf_file = os.path.join(os.path.dirname(__file__), "intron_compatibility.gtf")
        for locus_transcripts in parse_gtf(open(gtf_file)):
            isoform_graph = IsoformGraph.from_transcripts(locus_transcripts)
            isoform_graph.collapse(overhang_threshold=0)
            #write_dot(isoform_graph, "intron_compatibility.dot")
            self.compare_dot(isoform_graph, "intron_compatibility.dot")

    def testIncompatibleMerge(self):
        """test merge involving intron-incompatible paths"""
        gtf_file = os.path.join(os.path.dirname(__file__), "merge_incompatible.gtf")
        for locus_transcripts in parse_gtf(open(gtf_file)):
            isoform_graph = IsoformGraph.from_transcripts(locus_transcripts)
            #write_dot(isoform_graph, "merge_incompatible_b4.dot")
            isoform_graph.collapse()
            #write_dot(isoform_graph, "merge_incompatible_after.dot")
            #write_dot(isoform_graph, "merge_incompatible.dot")
            self.compare_dot(isoform_graph, "merge_incompatible.dot")

    def testMergeIdentical(self):
        gtf_file = os.path.join(os.path.dirname(__file__), "merge_identical.gtf")
        for locus_transcripts in parse_gtf(open(gtf_file)):
            isoform_graph = IsoformGraph.from_transcripts(locus_transcripts)
            isoform_graph.collapse()
            #write_dot(isoform_graph, "merge_identical.dot")
            self.compare_dot(isoform_graph, "merge_identical.dot")

    def testMergeSeparateStrands(self):
        gtf_file = os.path.join(os.path.dirname(__file__), "merge_separate_strands.gtf")
        for locus_transcripts in parse_gtf(open(gtf_file)):
            isoform_graph = IsoformGraph.from_transcripts(locus_transcripts)
            isoform_graph.collapse()
            #write_dot(isoform_graph, "merge_separate_strands.dot")
            self.compare_dot(isoform_graph, "merge_separate_strands.dot")

    def testMergeSeparateChromosomes(self):
        gtf_file = os.path.join(os.path.dirname(__file__), "merge_separate_chroms.gtf")
        for locus_transcripts in parse_gtf(open(gtf_file)):
            isoform_graph = IsoformGraph.from_transcripts(locus_transcripts) 
            isoform_graph.collapse()
            self.assertTrue(isoform_graph.G.number_of_nodes() == 1)

    def testMergeWithStrandlessExons(self):
        """test merging two stranded transcripts with an unstranded transcript"""
        gtf_file = os.path.join(os.path.dirname(__file__), "merge_nostrand.gtf")
        for locus_transcripts in parse_gtf(open(gtf_file)):
            isoform_graph = IsoformGraph.from_transcripts(locus_transcripts)            
            #write_dot(isoform_graph, "merge_nostrand_b4.dot")
            isoform_graph.collapse()
            #write_dot(isoform_graph, "merge_nostrand_after.dot")
            #write_dot(isoform_graph, "merge_nostrand.dot")
            self.compare_dot(isoform_graph, "merge_nostrand.dot")

    def testMultipleMerges(self):
        """test merge involving multiple paths"""
        gtf_file = os.path.join(os.path.dirname(__file__), "merge_multi.gtf")
        for locus_transcripts in parse_gtf(open(gtf_file)):
            isoform_graph = IsoformGraph.from_transcripts(locus_transcripts)        
            isoform_graph.collapse()
            #write_dot(isoform_graph, "merge_multi.dot")
            self.compare_dot(isoform_graph, "merge_multi.dot")

    def testMergeStrand(self):
        """test merge where only one exon has strand information"""
        gtf_file = os.path.join(os.path.dirname(__file__), "merge_strand.gtf")
        for locus_transcripts in parse_gtf(open(gtf_file)):
            isoform_graph = IsoformGraph.from_transcripts(locus_transcripts) 
            isoform_graph.collapse()
            #write_dot(isoform_graph, "merge_strand.dot")
            self.compare_dot(isoform_graph, "merge_strand.dot")

    def testMergeEncompassingPaths(self):
        """test merging transcripts where one transcript encompasses the other"""
        gtf_file = os.path.join(os.path.dirname(__file__), "merge_encompassing_path.gtf")
        for locus_transcripts in parse_gtf(open(gtf_file)):
            isoform_graph = IsoformGraph.from_transcripts(locus_transcripts)
            #write_dot(isoform_graph, "merge_encompassing_path_b4.dot")            
            isoform_graph.collapse()
            #write_dot(isoform_graph, "merge_encompassing_path_after.dot")            
            #write_dot(isoform_graph, "merge_encompassing_path.dot")      
            self.compare_dot(isoform_graph, "merge_encompassing_path.dot")

    def testMergePath(self):
        gtf_file = os.path.join(os.path.dirname(__file__), "merge_path.gtf")
        for locus_transcripts in parse_gtf(open(gtf_file)):
            isoform_graph = IsoformGraph.from_transcripts(locus_transcripts) 
            isoform_graph.collapse()
            #write_dot(isoform_graph, "merge_path.dot")
            self.compare_dot(isoform_graph, "merge_path.dot")

    def testMergePath2(self):
        """test a 3-way merge where a bridging exon joins two non-overlapping transcripts"""
        gtf_file = os.path.join(os.path.dirname(__file__), "merge_path2.gtf")
        for locus_transcripts in parse_gtf(open(gtf_file)):
            isoform_graph = IsoformGraph.from_transcripts(locus_transcripts)
            isoform_graph.collapse()
            #write_dot(isoform_graph, "merge_path2.dot")
            self.compare_dot(isoform_graph, "merge_path2.dot")
            
    def testMergePath3(self):
        """test merging exons that are touching but not overlapping"""
        gtf_file = os.path.join(os.path.dirname(__file__), "merge_path3.gtf")
        for locus_transcripts in parse_gtf(open(gtf_file)):
            isoform_graph = IsoformGraph.from_transcripts(locus_transcripts)
            #write_dot(isoform_graph, "merge_path3_b4.dot")
            isoform_graph.collapse()
            #write_dot(isoform_graph, "merge_path3_after.dot")
            #write_dot(isoform_graph, "merge_path3.dot")
            self.compare_dot(isoform_graph, "merge_path3.dot")

    def testMergePath4(self):
        """test merging two transcripts, where one is a truncated version of the other"""
        gtf_file = os.path.join(os.path.dirname(__file__), "merge_path4.gtf")
        for locus_transcripts in parse_gtf(open(gtf_file)):
            isoform_graph = IsoformGraph.from_transcripts(locus_transcripts)
            isoform_graph.collapse()
            #write_dot(isoform_graph, "merge_path4.dot")
            self.compare_dot(isoform_graph, "merge_path4.dot")

    def testMergePath5(self):
        """test merging two transcripts, where one is a truncated version of the other"""
        gtf_file = os.path.join(os.path.dirname(__file__), "merge_path5.gtf")
        for locus_transcripts in parse_gtf(open(gtf_file)):
            isoform_graph = IsoformGraph.from_transcripts(locus_transcripts)
            write_dot(isoform_graph, "merge_path5_b4.dot")
            isoform_graph.collapse()
            write_dot(isoform_graph, "merge_path5_after.dot")
            #write_dot(isoform_graph, "merge_path5.dot")
            nx.spring_layout(isoform_graph.G)
            nx.draw(isoform_graph.G)
            plt.show()            

            self.compare_dot(isoform_graph, "merge_path5.dot")

    def testMergePath6(self):
        """test a complex locus on chromosome 1 that gave errors due to problems with threading"""
        gtf_file = os.path.join(os.path.dirname(__file__), "merge_path6.gtf")
        for locus_transcripts in parse_gtf(open(gtf_file)):
            isoform_graph = IsoformGraph.from_transcripts(locus_transcripts)
            isoform_graph.collapse(trim=True, overhang_threshold=15)
            #write_dot(isoform_graph, "merge_path6_after.dot")
            #write_dot(isoform_graph, "merge_path6.dot")
            self.compare_dot(isoform_graph, "merge_path6.dot")
            break

    def testMergePath7(self):
        """test a complex locus on chromosome 1 that gave errors due to problems with threading"""
        gtf_file = os.path.join(os.path.dirname(__file__), "merge_path7.gtf")
        for locus_transcripts in parse_gtf(open(gtf_file)):
            isoform_graph = IsoformGraph.from_transcripts(locus_transcripts)
            isoform_graph.collapse(trim=True, overhang_threshold=15)
            #write_dot(isoform_graph, "merge_path7.dot")
            self.compare_dot(isoform_graph, "merge_path7.dot")

    def testMergePath8(self):
        """test a complex locus on chromosome 10 that caused problems"""
        gtf_file = os.path.join(os.path.dirname(__file__), "merge_path8.gtf")
        for locus_transcripts in parse_gtf(open(gtf_file)):
            isoform_graph = IsoformGraph.from_transcripts(locus_transcripts)
            #write_dot(isoform_graph, "merge_path8_b4.dot")            
            isoform_graph.collapse(trim=True, overhang_threshold=15)
            #write_dot(isoform_graph, "merge_path8_after.dot")            
            #write_dot(isoform_graph, "merge_path8.dot")            
            self.compare_dot(isoform_graph, "merge_path8.dot")

    def testMergePath9(self):
        """test a locus that causes an infinite merge loop"""
        gtf_file = os.path.join(os.path.dirname(__file__), "merge_path9.gtf")
        for locus_transcripts in parse_gtf(open(gtf_file)):
            isoform_graph = IsoformGraph.from_transcripts(locus_transcripts)
            isoform_graph.collapse(trim=True, overhang_threshold=15)
            #write_dot(isoform_graph, "merge_path9.dot")            
            self.compare_dot(isoform_graph, "merge_path9.dot")

    def testMergePath10(self):
        """test a locus that causes an infinite merge loop"""
        gtf_file = os.path.join(os.path.dirname(__file__), "merge_path10.gtf")
        for locus_transcripts in parse_gtf(open(gtf_file)):
            isoform_graph = IsoformGraph.from_transcripts(locus_transcripts)
            isoform_graph.collapse(trim=True, overhang_threshold=15)
            #write_dot(isoform_graph, "merge_path10.dot")            
            self.compare_dot(isoform_graph, "merge_path10.dot")

    def testTrim(self):
        """test trimming intron-incompatible paths"""
        gtf_file = os.path.join(os.path.dirname(__file__), "trim.gtf")
        for locus_transcripts in parse_gtf(open(gtf_file)):
            isoform_graph = IsoformGraph.from_transcripts(locus_transcripts)
            isoform_graph.collapse(trim=True, overhang_threshold=15)
            #write_dot(isoform_graph, "trim.dot")
            self.compare_dot(isoform_graph, "trim.dot")

    def testTrim2(self):
        """test a trimming example where a node can either be merged directly
        or trimmed first before merging"""
        gtf_file = os.path.join(os.path.dirname(__file__), "trim2.gtf")
        for locus_transcripts in parse_gtf(open(gtf_file)):
            isoform_graph = IsoformGraph.from_transcripts(locus_transcripts)
            #write_dot(isoform_graph, "trim2_b4.dot")            
            isoform_graph.collapse(trim=True, overhang_threshold=15)
            #write_dot(isoform_graph, "trim2_after.dot")            
            #write_dot(isoform_graph, "trim2.dot")
            self.compare_dot(isoform_graph, "trim2.dot")

    def testTrim3(self):
        """test locus that causes infinite merge loop"""
        gtf_file = os.path.join(os.path.dirname(__file__), "trim3.gtf")
        for locus_transcripts in parse_gtf(open(gtf_file)):
            isoform_graph = IsoformGraph.from_transcripts(locus_transcripts)
            #write_dot(isoform_graph, "trim3_b4.dot")            
            isoform_graph.collapse(trim=True, overhang_threshold=16)
            #write_dot(isoform_graph, "trim3_after.dot")            
            #write_dot(isoform_graph, "trim3.dot")            
            self.compare_dot(isoform_graph, "trim3.dot")
 
    def testThreading1(self):
        """ensure that leaf nodes that are 'absorbed' into the graph
        have their predecessor/successor nodes merged"""
        gtf_file = os.path.join(os.path.dirname(__file__), "threading1.gtf")
        for locus_transcripts in parse_gtf(open(gtf_file)):
            isoform_graph = IsoformGraph.from_transcripts(locus_transcripts)
            isoform_graph.collapse(trim=True, overhang_threshold=15)
            #write_dot(isoform_graph, "threading1.dot")  
            self.compare_dot(isoform_graph, "threading1.dot")

    def testClusterSingleExonTranscripts(self):
        gtf_file = os.path.join(os.path.dirname(__file__), "single_exon_transcripts1.gtf")
        for locus_transcripts in parse_gtf(open(gtf_file)):
            isoform_graph = IsoformGraph.from_transcripts(locus_transcripts)
            write_dot(isoform_graph, "single_exon_transcripts_b4.dot")  
            isoform_graph.collapse(trim=True, overhang_threshold=15)
            write_dot(isoform_graph, "single_exon_transcripts_after.dot")  
            #write_dot(isoform_graph, "single_exon_transcripts.dot")  
            self.compare_dot(isoform_graph, "single_exon_transcripts.dot")

    def testJoin1(self):
        """test whether compatible but unconnected exons can be linked to introns"""
        gtf_file = os.path.join(os.path.dirname(__file__), "join1.gtf")
        for locus_transcripts in parse_gtf(open(gtf_file)):
            isoform_graph = IsoformGraph.from_transcripts(locus_transcripts)
            isoform_graph.collapse(trim=True, overhang_threshold=15)
            #write_dot(isoform_graph, "join1.dot")  
            self.compare_dot(isoform_graph, "join1.dot")

    def testJoin2(self):
        """test join with unstranded node used by both + and - strand"""
        gtf_file = os.path.join(os.path.dirname(__file__), "join2.gtf")
        for locus_transcripts in parse_gtf(open(gtf_file)):
            isoform_graph = IsoformGraph.from_transcripts(locus_transcripts)
            isoform_graph.collapse(trim=True, overhang_threshold=15)
            #write_dot(isoform_graph, "join2.dot")  
            self.compare_dot(isoform_graph, "join2.dot")

    def testJoin3(self):
        """test join with unstranded node used by both + and - strand"""
        gtf_file = os.path.join(os.path.dirname(__file__), "join3.gtf")
        for locus_transcripts in parse_gtf(open(gtf_file)):
            isoform_graph = IsoformGraph.from_transcripts(locus_transcripts)
            #write_dot(isoform_graph, "join3_b4.dot")  
            isoform_graph.collapse(trim=True, overhang_threshold=15)
            #write_dot(isoform_graph, "join3_after.dot")  
            #write_dot(isoform_graph, "join3.dot")  
            self.compare_dot(isoform_graph, "join3.dot")



#            nx.spring_layout(isoform_graph.G)
#            nx.draw(isoform_graph.G)
#            plt.show()
#            nx.spring_layout(isoform_graph.G)
#            nx.draw(isoform_graph.G)
#            plt.show()
            #nx.spring_layout(isoform_graph.G)
            #nx.draw(isoform_graph.G)
            #plt.show()

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