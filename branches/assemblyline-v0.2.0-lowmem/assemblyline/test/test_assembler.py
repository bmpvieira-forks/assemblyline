'''
Created on Feb 20, 2013

@author: mkiyer
'''
import unittest

from assemblyline.lib.transcript import POS_STRAND, NEG_STRAND, Exon
from assemblyline.lib.assemble.transcript_graph import prune_transcript_graph
from assemblyline.lib.assemble.assembler import assemble_transcript_graph

from test_base import read_first_locus, get_transcript_graphs

class TestAssembler(unittest.TestCase):

    def test_assembler1(self):
        # setup correct transcripts
        PATH_ABCDE = tuple([Exon(0,100), Exon(200,300), Exon(400,500),Exon(600,700), Exon(800,900)])
        PATH_ACE = tuple([Exon(0,100), Exon(400,500), Exon(800,900)])
        PATH_ABCE = tuple([Exon(0,100), Exon(200,300), Exon(400,500), Exon(800,900)])
        PATH_ACDE = tuple([Exon(0,100), Exon(400,500),Exon(600,700), Exon(800,900)])
        # read transcripts
        transcripts = read_first_locus("assemble1.gtf", score_attr="score")
        GG = get_transcript_graphs(transcripts)
        G,tmap = GG[POS_STRAND]
        # set transcript scores
        tmap["ABCDE"].score = 2.0
        tmap["ACE"].score = 1.0
        tmap["ABCE"].score = 1.0
        tmap["ACDE"].score = 1.0
        # set assembly parameter
        kmax = 2
        # assemble
        GS = list(prune_transcript_graph(G, POS_STRAND, tmap,
                                         min_trim_length=0,
                                         trim_utr_fraction=0,
                                         trim_intron_fraction=0))
        Gsub, strand, partial_paths = GS[0]
        results = list(assemble_transcript_graph(Gsub, strand, partial_paths,
                                                 user_kmax=kmax,
                                                 ksensitivity=0,
                                                 fraction_major_path=0,
                                                 max_paths=1000))
        self.assertEquals(len(results), 2)
        self.assertEqual(tuple(results[0].path), PATH_ABCDE) 
        self.assertAlmostEqual(results[0].score, 3.0, places=3) 
        self.assertEqual(tuple(results[1].path), PATH_ACE) 
        self.assertAlmostEqual(results[1].score, 2.0, places=3) 
        # set transcript scores
        tmap["ABCDE"].score = 4.0
        tmap["ACE"].score = 3.0
        tmap["ABCE"].score = 2.0
        tmap["ACDE"].score = 1.0
        # set assembly parameter
        kmax = 3
        # assemble
        GS = list(prune_transcript_graph(G, POS_STRAND, tmap,
                                         min_trim_length=0,
                                         trim_utr_fraction=0,
                                         trim_intron_fraction=0))
        Gsub, strand, partial_paths = GS[0]
        results = list(assemble_transcript_graph(Gsub, strand, partial_paths,
                                                 user_kmax=kmax,
                                                 ksensitivity=0,
                                                 fraction_major_path=0,
                                                 max_paths=1000))
        self.assertEquals(len(results), 4)
        self.assertEqual(tuple(results[0].path), PATH_ABCDE) 
        self.assertAlmostEqual(results[0].score, 4.0, places=3) 
        self.assertEqual(tuple(results[1].path), PATH_ACE) 
        self.assertAlmostEqual(results[1].score, 3.0, places=3) 
        self.assertEqual(tuple(results[2].path), PATH_ABCE) 
        self.assertAlmostEqual(results[2].score, 2.0, places=3) 
        self.assertEqual(tuple(results[3].path), PATH_ACDE) 
        self.assertAlmostEqual(results[3].score, 1.0, places=3) 


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()