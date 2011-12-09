'''
Created on Feb 14, 2011

@author: mkiyer
'''
import logging
import unittest
import os
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from StringIO import StringIO

from assemblyline.lib.transcript import Exon, POS_STRAND, NEG_STRAND, NO_STRAND
from assemblyline.lib.transcript_graph import trim_left, trim_right, trim_transcript
from assemblyline.lib import transcript_graph
from assemblyline.lib.collapse_chains import collapse_contiguous_nodes, collapse_strand_specific_graph
from assemblyline.lib.transcript_parser import parse_gtf, cufflinks_attr_defs

from test_base import compare_dot, write_dot, convert_attrs_to_strings, read_gtf, make_transcript, get_gtf_path

logging.basicConfig(level=logging.DEBUG,
                    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")

class TestTranscriptGraph(unittest.TestCase):

    def test_basic_graph(self):
        """test that basic transcript parser works and builds graph"""        
        # plus strand
        test_basename = "basic_graph1"
        gtf_file = test_basename + ".gtf"
        dot_file = test_basename + ".dot"
        G = read_gtf(gtf_file)
        #write_dot(G, dot_file)
        self.assertTrue(compare_dot(G, dot_file))
        # minus strand
        test_basename = "basic_graph2"
        gtf_file = test_basename + ".gtf"
        dot_file = test_basename + ".dot"
        G = read_gtf(gtf_file)
        #write_dot(G, dot_file)
        self.assertTrue(compare_dot(G, dot_file))
        # unstranded multi-exon
        test_basename = "basic_graph3"
        gtf_file = test_basename + ".gtf"
        dot_file = test_basename + ".dot"
        G = read_gtf(gtf_file)
        #write_dot(G, dot_file)
        self.assertTrue(compare_dot(G, dot_file))
        #nx.spring_layout(G)
        #nx.draw(G)
        #plt.show()

    def test_basic_split(self):
        """test that initial graph construction correctly splits exons"""
        # plus strand
        test_basename = "basic_split1"
        gtf_file = test_basename + ".gtf"
        dot_file = test_basename + ".dot"
        G = read_gtf(gtf_file)
        #write_dot(G, dot_file)
        self.assertTrue(compare_dot(G, dot_file))       
        # minus strand
        test_basename = "basic_split2"
        gtf_file = test_basename + ".gtf"
        dot_file = test_basename + ".dot"
        G = read_gtf(gtf_file)
        #write_dot(G, dot_file)
        self.assertTrue(compare_dot(G, dot_file))
        # unstranded
        test_basename = "basic_split3"
        gtf_file = test_basename + ".gtf"
        dot_file = test_basename + ".dot"
        G = read_gtf(gtf_file)
        #write_dot(G, dot_file)
        self.assertTrue(compare_dot(G, dot_file))
        #nx.spring_layout(G)
        #nx.draw(G)
        #plt.show()

    def test_complex_split(self):
        """test some more complicated exon splitting cases"""
        # only one strand
        test_basename = "complex_split1"
        gtf_file = test_basename + ".gtf"
        dot_file = test_basename + ".dot"
        G = read_gtf(gtf_file)
        #write_dot(G, dot_file)
        self.assertTrue(compare_dot(G, dot_file))
        # duplicate plus and minus strand overlapping
        test_basename = "complex_split2"
        gtf_file = test_basename + ".gtf"
        dot_file = test_basename + ".dot"
        G = read_gtf(gtf_file)
        #write_dot(G, dot_file)
        self.assertTrue(compare_dot(G, dot_file))
        # add in unstranded transcripts
        test_basename = "complex_split3"
        gtf_file = test_basename + ".gtf"
        dot_file = test_basename + ".dot"
        G = read_gtf(gtf_file)
        #write_dot(G, dot_file)
        self.assertTrue(compare_dot(G, dot_file))
        #nx.spring_layout(G)
        #nx.draw(G)
        #plt.show()

    def test_antisense(self):
        """test some examples of antisense transcripts"""
        # exons do not overlap
        test_basename = "antisense1"
        gtf_file = test_basename + ".gtf"
        dot_file = test_basename + ".dot"
        G = read_gtf(gtf_file)
        #write_dot(G, dot_file)
        self.assertTrue(compare_dot(G, dot_file))
        # exons partially overlap
        test_basename = "antisense2"
        gtf_file = test_basename + ".gtf"
        dot_file = test_basename + ".dot"
        G = read_gtf(gtf_file)
        #write_dot(G, dot_file)
        self.assertTrue(compare_dot(G, dot_file))
        # one exon shared among strands
        test_basename = "antisense3"
        gtf_file = test_basename + ".gtf"
        dot_file = test_basename + ".dot"
        G = read_gtf(gtf_file)
        #write_dot(G, dot_file)
        self.assertTrue(compare_dot(G, dot_file))
        #nx.spring_layout(G)
        #nx.draw(G)
        #plt.show()
       
#    def test_example1(self):
#        """test some examples from actual data"""
#        # exons do not overlap
#        test_basename = "ex1"
#        gtf_file = test_basename + ".gtf"
#        dot_file = test_basename + ".dot"
#        txgraph = read_gtf(gtf_file)
#        write_dot(txgraph, dot_file)
#        self.assertTrue(compare_dot(txgraph, dot_file))        


class TestTrimming(unittest.TestCase):

    def test_trim_left_right(self):
        # make up some intron start/end intervals
        intron_starts = [100, 400, 700]
        intron_ends = [200, 500, 800]
        # try trimming fake exon left
        overhang_threshold = 10
        x = trim_left(195, 220, intron_starts, intron_ends, overhang_threshold)
        self.assertEqual(x, 200)
        x = trim_left(191, 220, intron_starts, intron_ends, overhang_threshold)
        self.assertEqual(x, 200)
        x = trim_left(190, 220, intron_starts, intron_ends, overhang_threshold)
        self.assertEqual(x, 200)
        x = trim_left(189, 220, intron_starts, intron_ends, overhang_threshold)
        self.assertEqual(x, 189)
        # try trimming fake exon right
        x = trim_right(300, 405, intron_starts, intron_ends, overhang_threshold)
        self.assertEqual(x, 400)
        x = trim_right(300, 409, intron_starts, intron_ends, overhang_threshold)
        self.assertEqual(x, 400)
        x = trim_right(300, 410, intron_starts, intron_ends, overhang_threshold)
        self.assertEqual(x, 400)
        x = trim_right(300, 411, intron_starts, intron_ends, overhang_threshold)
        self.assertEqual(x, 411)
        
    def test_trim_transcript(self):        
        # make up some intron start/end intervals
        intron_starts = [100, 400, 700]
        intron_ends = [200, 500, 800]
        overhang_threshold = 10        
        # trim some fake transcripts
        s,e = trim_transcript([Exon(195,405)], intron_starts, intron_ends, 
                              overhang_threshold)
        self.assertEquals((s,e), (200,400))
        s,e = trim_transcript([Exon(190,410)], intron_starts, intron_ends, 
                              overhang_threshold)
        self.assertEquals((s,e), (200,400))
        s,e = trim_transcript([Exon(189,411)], intron_starts, intron_ends, 
                              overhang_threshold)
        self.assertEquals((s,e), (189,411))
        # create situation where transcript trimming conflicts
        intron_starts = [100, 190, 400]
        intron_ends = [200, 500]
        overhang_threshold = 20
        s,e = trim_transcript([Exon(188,201)], intron_starts, intron_ends, 
                              overhang_threshold)
        self.assertEquals((s,e), (188,190))
        # behavior in the case of a tie is to trim right (arbitrary)
        s,e = trim_transcript([Exon(189,201)], intron_starts, intron_ends, 
                              overhang_threshold)
        self.assertEquals((s,e), (200,201))
        s,e = trim_transcript([Exon(189,202)], intron_starts, intron_ends, 
                              overhang_threshold)
        self.assertEquals((s,e), (200,202))


class TestRedistributeDensity(unittest.TestCase):

    def load_gtf(self, gtf_file):
        for locus_transcripts in parse_gtf(open(get_gtf_path(gtf_file)), cufflinks_attr_defs):
            G = transcript_graph.add_transcripts_to_graph(locus_transcripts)
            break
        return G, locus_transcripts

    def test_redistribute(self):
        """
        ensure that coverage from unstranded transcripts is allocated
        proportionally to the coverage from overlapping stranded 
        transcripts
        """
        #
        # case where an unstranded transcript lies in UTR region and
        # overlaps both a + and - strand transcript
        #
        gtf_file = "redistribute1.gtf"
        G, transcripts = self.load_gtf(gtf_file)
        # check density arrays at key nodes before redistribution
        a = G.node[Exon(950,1000)]['strand_density']
        self.assertTrue(np.array_equal(a, np.array((1.0, 0.0, 1.0))))
        a = G.node[Exon(1000,1100)]['strand_density']
        self.assertTrue(np.array_equal(a, np.array((0.0, 0.0, 1.0))))
        a = G.node[Exon(1100,1150)]['strand_density']
        self.assertTrue(np.array_equal(a, np.array((0.0, 1.0, 1.0))))
        # redistribute coverage
        transcript_graph.redistribute_unstranded_density(G, transcripts)
        # check density after redistribution
        a = G.node[Exon(950,1000)]['strand_density']
        self.assertTrue(np.array_equal(a, np.array((1.5, 0.5, 0.0))))
        a = G.node[Exon(1000,1100)]['strand_density']
        self.assertTrue(np.array_equal(a, np.array((0.5, 0.5, 0.0))))
        a = G.node[Exon(1100,1150)]['strand_density']
        self.assertTrue(np.array_equal(a, np.array((0.5, 1.5, 0.0))))
        #
        # case where several unstranded transcripts are alone and do not 
        # overlap stranded transcripts
        #
        gtf_file = "redistribute2.gtf"
        G, transcripts = self.load_gtf(gtf_file)
        # check density arrays before redistribution
        for n,d in G.nodes_iter(data=True):
            self.assertAlmostEqual(d['strand_density'][POS_STRAND], 0)
            self.assertAlmostEqual(d['strand_density'][NEG_STRAND], 0)
            self.assertTrue(d['strand_density'][NO_STRAND] > 0)
        # redistribute coverage
        transcript_graph.redistribute_unstranded_density(G, transcripts)
        # check density arrays after redistribution
        for n,d in G.nodes_iter(data=True):
            self.assertAlmostEqual(d['strand_density'][POS_STRAND], 0)
            self.assertAlmostEqual(d['strand_density'][NEG_STRAND], 0)
            self.assertTrue(d['strand_density'][NO_STRAND] > 0)        
        #
        # case with a chain of unstranded transcripts and only one of
        # them overlaps a stranded transcript 
        #
        gtf_file = "redistribute3.gtf"
        G, transcripts = self.load_gtf(gtf_file)
        # check density before redistribution
        exons_overlapping_pos = [Exon(200,210),Exon(290,300)]
        exons_not_overlapping = [Exon(50,110), Exon(110,150), Exon(150,200), 
                                 Exon(300,350), Exon(350,390), Exon(390,450)]
        for e in exons_overlapping_pos:        
            a = G.node[e]['strand_density']
            self.assertTrue(np.array_equal(a, np.array((1.0, 0.0, 1.0))))
        for e in exons_not_overlapping: 
            a = G.node[e]['strand_density']
            self.assertTrue(a[POS_STRAND] == 0)
            self.assertTrue(a[NO_STRAND] > 0.0)
        # redistribute
        transcript_graph.redistribute_unstranded_density(G, transcripts)
        GG = transcript_graph.create_strand_specific_graphs(G)
        G = collapse_strand_specific_graph(GG[POS_STRAND])
        # density should be correct
        self.assertTrue(G.number_of_nodes() == 1)
        self.assertTrue(G.node[Exon(50,450)]['density'], 500./400.)


class TestCollapseChains(unittest.TestCase):
    
    def load_gtf(self, gtf_file):
        for locus_transcripts in parse_gtf(open(get_gtf_path(gtf_file)), cufflinks_attr_defs):
            G = transcript_graph.add_transcripts_to_graph(locus_transcripts)
            break
        # redistribute coverage
        transcript_graph.redistribute_unstranded_density(G, locus_transcripts)
        # partition into single-stranded graphs
        GG = transcript_graph.create_strand_specific_graphs(G)
        return GG
        
    def test_collapse_chains(self):
        # test on purely unstranded data
        gtf_file = "redistribute2.gtf"
        GG = self.load_gtf(gtf_file)
        G = GG[NO_STRAND]
        self.assertEqual(G.number_of_nodes(), 3)
        H = collapse_strand_specific_graph(GG[NO_STRAND])
        self.assertEqual(H.number_of_nodes(), 1)
        self.assertAlmostEqual(H.node[Exon(0,1000)]['density'], 1.6)
        # test on single exon (+) strand transcript with several
        # overlapping unknown-strand transcripts and a (-) strand
        # transcript at end that should not interfere
        gtf_file = "collapse1.gtf"
        GG = self.load_gtf(gtf_file)
        Hfwd = collapse_strand_specific_graph(GG[POS_STRAND])
        Hrev = collapse_strand_specific_graph(GG[NEG_STRAND])
        Hunknown = collapse_strand_specific_graph(GG[NO_STRAND])
        self.assertEqual(Hfwd.number_of_nodes(), 1)
        self.assertAlmostEqual(Hfwd.node[Exon(0,1000)]['density'], 1.5)
        self.assertEqual(Hrev.number_of_nodes(), 1)
        self.assertAlmostEqual(Hrev.node[Exon(900,1000)]['density'], 1.0)
        self.assertEqual(Hunknown.number_of_nodes(), 0)


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()