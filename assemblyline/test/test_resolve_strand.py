'''
Created on Feb 20, 2013

@author: mkiyer
'''
import unittest

from assemblyline.lib.assemble.transcript_graph import \
    partition_transcripts_by_strand
from assemblyline.lib.transcript import POS_STRAND, NEG_STRAND, NO_STRAND
from test_base import read_first_locus

class TestResolveStrand(unittest.TestCase):
    def test_resolve_strand_ref(self):
        transcripts = read_first_locus("resolve_strand_ref1.gtf", score_attr="FPKM")
        tdict = dict((t.attrs['transcript_id'],t) for t in transcripts)
        partition_transcripts_by_strand(transcripts)
        # check resolved strands
        t = tdict["T1"]
        self.assertTrue(t.strand == POS_STRAND)
        t = tdict["T2"]
        self.assertTrue(t.strand == NEG_STRAND)
        # equal overlap on both strands, default to positive
        t = tdict["T3"]
        self.assertTrue(t.strand == POS_STRAND)
        # more positive strand overlap
        t = tdict["T4"]
        self.assertTrue(t.strand == POS_STRAND)
        # more negative strand overlap
        t = tdict["T5"]
        self.assertTrue(t.strand == NEG_STRAND)
        return
    
    def test_resolve_strand_scores1(self):
        transcripts = read_first_locus("resolve_strand_scores1.gtf", score_attr="FPKM")
        tdict = dict((t.attrs['transcript_id'],t) for t in transcripts)
        partition_transcripts_by_strand(transcripts)
        # check strands
        t = tdict["AA"]
        self.assertTrue(t.strand == NO_STRAND)
        t = tdict["BB"]
        self.assertTrue(t.strand == NO_STRAND)

    def test_resolve_strand_scores2(self):
        transcripts = read_first_locus("resolve_strand_scores2.gtf", score_attr="FPKM")
        tdict = dict((t.attrs['transcript_id'],t) for t in transcripts)
        partition_transcripts_by_strand(transcripts)
        # check strands
        t = tdict["AA"]
        self.assertTrue(t.strand == POS_STRAND)
        t = tdict["BB"]
        self.assertTrue(t.strand == POS_STRAND)
        t = tdict["CC"]
        self.assertTrue(t.strand == NEG_STRAND)
        t = tdict["DD"]
        self.assertTrue(t.strand == NEG_STRAND)

    def test_resolve_strand_scores3(self):
        transcripts = read_first_locus("resolve_strand_scores3.gtf", score_attr="FPKM")
        tdict = dict((t.attrs['transcript_id'],t) for t in transcripts)
        partition_transcripts_by_strand(transcripts)
        # check strands
        t = tdict['CC']
        self.assertTrue(t.strand == POS_STRAND)
        t = tdict['DD']
        self.assertTrue(t.strand == NEG_STRAND)
        t = tdict['EE']
        self.assertTrue(t.strand == POS_STRAND)
        # increase score on negative strand
        transcripts = read_first_locus("resolve_strand_scores3.gtf", score_attr="FPKM")
        tdict = dict((t.attrs['transcript_id'],t) for t in transcripts)
        tdict['DD'].score += 1.0
        partition_transcripts_by_strand(transcripts)
        # check strands
        t = tdict['EE']
        self.assertTrue(t.strand == NEG_STRAND)


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()