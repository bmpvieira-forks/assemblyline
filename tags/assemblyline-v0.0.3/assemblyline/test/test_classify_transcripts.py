'''
Created on Nov 30, 2011

@author: mkiyer
'''
import unittest

from assemblyline.lib.transcript import Transcript, Exon, strand_str_to_int
from assemblyline.pipeline.classify_transcripts import annotate_recurrence_and_density, RECUR, DENSITY

def make_transcript(chrom='chr1', strand='.', 
                    id='Alib.1', sample='A', cohort='AA', 
                    fpkm=0.0, cov=0.0, frac=0.0, exons=None,
                    attrs=None):
    start = min(e[0] for e in exons)
    end = max(e[1] for e in exons)
    t = Transcript()
    t.chrom = chrom
    t.start = start
    t.end = end
    t.strand = strand_str_to_int(strand)
    t.id = id
    t.library = id.split(".")[0]
    t.sample = sample
    t.cohort = cohort
    t.fpkm = fpkm
    t.cov = cov
    t.frac = frac
    t.exons = [Exon(e[0], e[1]) for e in exons]
    t.attrs = {}
    return t

t1 = make_transcript(sample="A", strand="+", frac=1.0, fpkm=1.0, cov=1.0, exons=[(0,2), (5,7), (9,10)])
t2 = make_transcript(sample="A", strand="+", frac=1.0, fpkm=1.0, cov=1.0, exons=[(0,2), (6,7), (9,10)])
t3 = make_transcript(sample="B", strand="+", frac=1.0, fpkm=1.0, cov=1.0, exons=[(0,7)])

t4 = make_transcript(sample="A", strand="-", frac=1.0, fpkm=1.0, cov=1.0, exons=[(0,2), (5,7), (9,10)])
t5 = make_transcript(sample="A", strand="-", frac=1.0, fpkm=1.0, cov=1.0, exons=[(0,2), (6,7), (9,10)])
t6 = make_transcript(sample="B", strand="-", frac=1.0, fpkm=1.0, cov=1.0, exons=[(0,7)])

class TestAnnotateRecurrenceDensity(unittest.TestCase):
    
    def test_stranded_transcripts(self):
        # test just one sample
        libs_per_sample = {'A': 1}
        annotate_recurrence_and_density([t1], libs_per_sample)
        self.assertAlmostEqual(t1.attrs[RECUR], 1.0)
        # test one sample with two transcripts
        libs_per_sample = {'A': 1}
        annotate_recurrence_and_density([t1, t2], libs_per_sample)
        self.assertAlmostEqual(t1.attrs[RECUR], 1.0)
        self.assertAlmostEqual(t2.attrs[RECUR], 1.0)
        # test one sample with two transcripts on opposite strands
        libs_per_sample = {'A': 1}
        annotate_recurrence_and_density([t1, t4], libs_per_sample)
        self.assertAlmostEqual(t1.attrs[RECUR], 1.0)
        self.assertAlmostEqual(t2.attrs[RECUR], 1.0)        
        # test two samples
        libs_per_sample = {'A': 1, 'B': 1}
        annotate_recurrence_and_density([t1,t2,t3], libs_per_sample)
        self.assertAlmostEqual(t1.attrs[RECUR], 1.80)
        self.assertAlmostEqual(t2.attrs[RECUR], 1.75)
        self.assertAlmostEqual(t3.attrs[RECUR], 1.0 + 4.0/7.0)
        # test two samples with transcripts on opposite strands
        libs_per_sample = {'A': 1, 'B': 1}
        annotate_recurrence_and_density([t1,t2,t3,t4,t5,t6], libs_per_sample)
        self.assertAlmostEqual(t1.attrs[RECUR], 1.80)
        self.assertAlmostEqual(t2.attrs[RECUR], 1.75)
        self.assertAlmostEqual(t3.attrs[RECUR], 1.0 + 4.0/7.0)
        self.assertAlmostEqual(t4.attrs[RECUR], 1.80)
        self.assertAlmostEqual(t5.attrs[RECUR], 1.75)
        self.assertAlmostEqual(t6.attrs[RECUR], 1.0 + 4.0/7.0)
        
    def test_unstranded_transcripts(self):
        pass


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()