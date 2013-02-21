'''
Created on Nov 8, 2012

@author: mkiyer
'''
import unittest

from assemblyline.lib.assemble.trim import trim_graph
from assemblyline.lib.transcript import POS_STRAND, NEG_STRAND, Exon

from test_base import read_first_locus, get_transcript_graphs

class TestTrim(unittest.TestCase):

    def test_trim_bidir(self):
        transcripts = read_first_locus("trim_bidir1.gtf", score_attr="FPKM")
        GG = get_transcript_graphs(transcripts)
        G,tmap = GG[POS_STRAND]
        # trim at three different thresholds
        trim_nodes = trim_graph(G, POS_STRAND,
                                min_trim_length=0, 
                                trim_utr_fraction=0.015,
                                trim_intron_fraction=0.0)
        correct = set([Exon(0,100), Exon(900,1000)])
        self.assertTrue(trim_nodes == correct)
        trim_nodes = trim_graph(G, POS_STRAND,
                                min_trim_length=0, 
                                trim_utr_fraction=0.11,
                                trim_intron_fraction=0.0)
        correct = set([Exon(0,100), Exon(900,1000), 
                       Exon(100,200), Exon(800,900)])
        self.assertTrue(trim_nodes == correct)
        trim_nodes = trim_graph(G, POS_STRAND,
                                min_trim_length=0, 
                                trim_utr_fraction=0.26,
                                trim_intron_fraction=0.0)
        correct = set([Exon(0,100), Exon(900,1000), 
                       Exon(100,200), Exon(800,900),
                       Exon(200,300), Exon(700,800)])
        self.assertTrue(trim_nodes == correct)
        # flip sign of transcripts and try again
        for t in transcripts:
            t.strand = NEG_STRAND
        GG = get_transcript_graphs(transcripts)
        G,tmap = GG[NEG_STRAND]        
        # trim at three different thresholds
        trim_nodes = trim_graph(G, NEG_STRAND,
                                min_trim_length=0, 
                                trim_utr_fraction=0.015,
                                trim_intron_fraction=0.0)
        correct = set([Exon(0,100), Exon(900,1000)])
        self.assertTrue(trim_nodes == correct)
        trim_nodes = trim_graph(G, NEG_STRAND,
                                min_trim_length=0, 
                                trim_utr_fraction=0.11,
                                trim_intron_fraction=0.0)
        correct = set([Exon(0,100), Exon(900,1000), 
                       Exon(100,200), Exon(800,900)])
        self.assertTrue(trim_nodes == correct)
        trim_nodes = trim_graph(G, NEG_STRAND,
                                min_trim_length=0, 
                                trim_utr_fraction=0.26,
                                trim_intron_fraction=0.0)
        correct = set([Exon(0,100), Exon(900,1000), 
                       Exon(100,200), Exon(800,900),
                       Exon(200,300), Exon(700,800)])
        self.assertTrue(trim_nodes == correct)

    def test_trim_intron_bidir(self):
        transcripts = read_first_locus("trim_intron_bidir1.gtf", score_attr="FPKM")
        GG = get_transcript_graphs(transcripts)
        G,tmap = GG[POS_STRAND]      
        # trim at different thresholds
        trim_nodes = trim_graph(G, POS_STRAND,
                                min_trim_length=0, 
                                trim_utr_fraction=0.0,
                                trim_intron_fraction=0.001)
        correct = set()
        self.assertTrue(trim_nodes == correct)
        trim_nodes = trim_graph(G, POS_STRAND,
                                min_trim_length=0, 
                                trim_utr_fraction=0.0,
                                trim_intron_fraction=0.025)
        correct = set([Exon(1900, 2000), Exon(1000, 1100)])
        self.assertTrue(trim_nodes == correct)
        trim_nodes = trim_graph(G, POS_STRAND,
                                min_trim_length=0, 
                                trim_utr_fraction=0.0,
                                trim_intron_fraction=0.2)
        correct = set([Exon(1900, 2000), 
                       Exon(1100, 1200),
                       Exon(1800, 1900), 
                       Exon(1000, 1100)])
        self.assertTrue(trim_nodes == correct)
        trim_nodes = trim_graph(G, POS_STRAND,
                                min_trim_length=0, 
                                trim_utr_fraction=0.0,
                                trim_intron_fraction=0.25)
        correct = set([Exon(1900, 2000), 
                       Exon(1100, 1200),
                       Exon(1200, 1300), 
                       Exon(1700, 1800), 
                       Exon(1800, 1900), 
                       Exon(1000, 1100)])
        self.assertTrue(trim_nodes == correct)
        # flip sign of transcripts and try again
        for t in transcripts:
            t.strand = NEG_STRAND
        GG = get_transcript_graphs(transcripts)
        G,tmap = GG[NEG_STRAND]
        # trim at different thresholds
        trim_nodes = trim_graph(G, NEG_STRAND,
                                min_trim_length=0, 
                                trim_utr_fraction=0.0,
                                trim_intron_fraction=0.001)
        correct = set()
        self.assertTrue(trim_nodes == correct)
        trim_nodes = trim_graph(G, NEG_STRAND,
                                min_trim_length=0, 
                                trim_utr_fraction=0.0,
                                trim_intron_fraction=0.025)
        correct = set([Exon(1900, 2000), Exon(1000, 1100)])
        self.assertTrue(trim_nodes == correct)
        trim_nodes = trim_graph(G, NEG_STRAND,
                                min_trim_length=0, 
                                trim_utr_fraction=0.0,
                                trim_intron_fraction=0.2)
        correct = set([Exon(1900, 2000), 
                       Exon(1100, 1200),
                       Exon(1800, 1900), 
                       Exon(1000, 1100)])
        self.assertTrue(trim_nodes == correct)
        trim_nodes = trim_graph(G, NEG_STRAND,
                                min_trim_length=0, 
                                trim_utr_fraction=0.0,
                                trim_intron_fraction=0.25)
        correct = set([Exon(1900, 2000), 
                       Exon(1100, 1200),
                       Exon(1200, 1300), 
                       Exon(1700, 1800), 
                       Exon(1800, 1900), 
                       Exon(1000, 1100)])
        self.assertTrue(trim_nodes == correct)


    def test_trim_utr(self):
        transcripts = read_first_locus("trim_utr1.gtf", score_attr="FPKM")
        GG = get_transcript_graphs(transcripts)
        G,tmap = GG[POS_STRAND]
        # trim at different thresholds
        trim_nodes = trim_graph(G, POS_STRAND,
                                min_trim_length=0, 
                                trim_utr_fraction=0.009,
                                trim_intron_fraction=0.0)
        correct = set()
        self.assertTrue(trim_nodes == correct)        
        trim_nodes = trim_graph(G, POS_STRAND,
                                min_trim_length=0, 
                                trim_utr_fraction=0.015,
                                trim_intron_fraction=0.0)
        correct = set([Exon(1000,1100)])
        self.assertTrue(trim_nodes == correct)
        trim_nodes = trim_graph(G, POS_STRAND,
                                min_trim_length=0, 
                                trim_utr_fraction=0.11,
                                trim_intron_fraction=0.0)
        correct = set([Exon(1000,1100), Exon(1100,1200)])
        self.assertTrue(trim_nodes == correct)
        trim_nodes = trim_graph(G, POS_STRAND,
                                min_trim_length=0, 
                                trim_utr_fraction=0.26,
                                trim_intron_fraction=0.0)
        correct = set([Exon(1000,1100), Exon(1100,1200),
                       Exon(1200,1300)])
        self.assertTrue(trim_nodes == correct)
        trim_nodes = trim_graph(G, POS_STRAND,
                                min_trim_length=0, 
                                trim_utr_fraction=1.0,
                                trim_intron_fraction=0.0)
        correct = set([Exon(1000,1100), Exon(1100,1200),
                       Exon(1200,1300)])
        self.assertTrue(trim_nodes == correct)        
        # flip sign of transcripts and try again
        for t in transcripts:
            t.strand = NEG_STRAND
        GG = get_transcript_graphs(transcripts)
        G,tmap = GG[NEG_STRAND]
        # trim at different thresholds
        trim_nodes = trim_graph(G, NEG_STRAND,
                                min_trim_length=0, 
                                trim_utr_fraction=0.009,
                                trim_intron_fraction=0.0)
        correct = set()
        self.assertTrue(trim_nodes == correct)  
        trim_nodes = trim_graph(G, NEG_STRAND,
                                min_trim_length=0, 
                                trim_utr_fraction=0.015,
                                trim_intron_fraction=0.0)
        correct = set([Exon(1000,1100)])
        self.assertTrue(trim_nodes == correct)
        trim_nodes = trim_graph(G, NEG_STRAND,
                                min_trim_length=0, 
                                trim_utr_fraction=0.11,
                                trim_intron_fraction=0.0)
        correct = set([Exon(1000,1100), Exon(1100,1200)])
        self.assertTrue(trim_nodes == correct)
        trim_nodes = trim_graph(G, NEG_STRAND,
                                min_trim_length=0, 
                                trim_utr_fraction=0.26,
                                trim_intron_fraction=0.0)
        correct = set([Exon(1000,1100), Exon(1100,1200),
                       Exon(1200,1300)])
        self.assertTrue(trim_nodes == correct)
        trim_nodes = trim_graph(G, NEG_STRAND,
                                min_trim_length=0, 
                                trim_utr_fraction=1.0,
                                trim_intron_fraction=0.0)
        correct = set([Exon(1000,1100), Exon(1100,1200),
                       Exon(1200,1300)])
        self.assertTrue(trim_nodes == correct)  

    def test_trim_intronic_utr(self):
        transcripts = read_first_locus("trim_intron_utr1.gtf", score_attr="FPKM")
        GG = get_transcript_graphs(transcripts)
        G,tmap = GG[POS_STRAND]       
        # trim at different thresholds
        trim_nodes = trim_graph(G, POS_STRAND,
                                min_trim_length=0, 
                                trim_utr_fraction=0.0,
                                trim_intron_fraction=0.001)
        correct = set()
        self.assertTrue(trim_nodes == correct)
        trim_nodes = trim_graph(G, POS_STRAND,
                                min_trim_length=0, 
                                trim_utr_fraction=0.0,
                                trim_intron_fraction=0.011)
        correct = set([Exon(1000,1100)])
        self.assertTrue(trim_nodes == correct)
        trim_nodes = trim_graph(G, POS_STRAND,
                                min_trim_length=0, 
                                trim_utr_fraction=0.0,
                                trim_intron_fraction=0.055)
        correct = set([Exon(1000,1100), Exon(1100,1200)])
        self.assertTrue(trim_nodes == correct)
        trim_nodes = trim_graph(G, POS_STRAND,
                                min_trim_length=0, 
                                trim_utr_fraction=0.0,
                                trim_intron_fraction=0.15)
        correct = set([Exon(1000,1100), Exon(1100,1200),
                       Exon(1200,1300)])
        self.assertTrue(trim_nodes == correct)
        trim_nodes = trim_graph(G, POS_STRAND,
                                min_trim_length=0, 
                                trim_utr_fraction=0.0,
                                trim_intron_fraction=1.0)
        correct = set([Exon(1000,1100), Exon(1100,1200),
                       Exon(1200,1300), Exon(1300,1500)])
        self.assertTrue(trim_nodes == correct)        
        # flip sign of transcripts and try again
        for t in transcripts:
            t.strand = NEG_STRAND
        GG = get_transcript_graphs(transcripts)
        G,tmap = GG[NEG_STRAND]
        # trim at different thresholds
        trim_nodes = trim_graph(G, NEG_STRAND,
                                min_trim_length=0, 
                                trim_utr_fraction=0.0,
                                trim_intron_fraction=0.001)
        correct = set()
        self.assertTrue(trim_nodes == correct)
        trim_nodes = trim_graph(G, NEG_STRAND,
                                min_trim_length=0, 
                                trim_utr_fraction=0.0,
                                trim_intron_fraction=0.011)
        correct = set([Exon(1000,1100)])
        self.assertTrue(trim_nodes == correct)
        trim_nodes = trim_graph(G, NEG_STRAND,
                                min_trim_length=0, 
                                trim_utr_fraction=0.0,
                                trim_intron_fraction=0.055)
        correct = set([Exon(1000,1100), Exon(1100,1200)])
        self.assertTrue(trim_nodes == correct)
        trim_nodes = trim_graph(G, NEG_STRAND,
                                min_trim_length=0, 
                                trim_utr_fraction=0.0,
                                trim_intron_fraction=0.15)
        correct = set([Exon(1000,1100), Exon(1100,1200),
                       Exon(1200,1300)])
        self.assertTrue(trim_nodes == correct)
        trim_nodes = trim_graph(G, NEG_STRAND,
                                min_trim_length=0, 
                                trim_utr_fraction=0.0,
                                trim_intron_fraction=1.0)
        correct = set([Exon(1000,1100), Exon(1100,1200),
                       Exon(1200,1300), Exon(1300,1500)])
        self.assertTrue(trim_nodes == correct)

    def test_trim_intron_retention(self):
        transcripts = read_first_locus("trim_intron_retention1.gtf", score_attr="FPKM")
        GG = get_transcript_graphs(transcripts)
        G,tmap = GG[POS_STRAND]
        # trim at different thresholds
        trim_nodes = trim_graph(G, POS_STRAND,
                                min_trim_length=0, 
                                trim_utr_fraction=0.0,
                                trim_intron_fraction=0.01)
        correct = set()
        self.assertTrue(trim_nodes == correct)
        trim_nodes = trim_graph(G, POS_STRAND,
                                min_trim_length=0, 
                                trim_utr_fraction=0.0,
                                trim_intron_fraction=0.11)
        correct = set([Exon(500,1500)])
        self.assertTrue(trim_nodes == correct)
        trim_nodes = trim_graph(G, POS_STRAND,
                                min_trim_length=0, 
                                trim_utr_fraction=0.0,
                                trim_intron_fraction=0.21)
        correct = set([Exon(500,1500), Exon(2000,9000)])
        self.assertTrue(trim_nodes == correct)
        trim_nodes = trim_graph(G, POS_STRAND,
                                min_trim_length=0, 
                                trim_utr_fraction=0.0,
                                trim_intron_fraction=1.0)
        correct = set([Exon(500,1500), Exon(2000,9000)])
        self.assertTrue(trim_nodes == correct)


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()