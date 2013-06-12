'''
Created on Feb 12, 2013

@author: mkiyer
'''
import unittest

# project imports
from assemblyline.pipeline.annotate_transcripts import annotate_locus
from assemblyline.lib.base import GTFAttr

# local imports
from test_base import read_first_locus

class TestAnnotate(unittest.TestCase):

    def test_categories(self):
        transcripts = read_first_locus("annotate_category1.gtf")
        t_dict = dict((t.attrs['transcript_id'],t) for t in transcripts)
        annotate_locus(transcripts, 
                       gtf_sample_attr="sample_id")
        # intronic same strand
        self.assertTrue(t_dict['T2'].attrs[GTFAttr.CATEGORY] == 2)
        # intronic opposite strand
        self.assertTrue(t_dict['T3'].attrs[GTFAttr.CATEGORY] == 3)
        # intronic ambiguous
        self.assertTrue(t_dict['T6'].attrs[GTFAttr.CATEGORY] == 4)
        # interleaving
        self.assertTrue(t_dict['T4'].attrs[GTFAttr.CATEGORY] == 5)
        # interleaving
        self.assertTrue(t_dict['T5'].attrs[GTFAttr.CATEGORY] == 3)
        # opp strand overlap (no introns)
        t = t_dict['T7']
        self.assertTrue(t.attrs[GTFAttr.CATEGORY] == 1)
        self.assertTrue(t.attrs[GTFAttr.ANN_REF_ID] == 'T1')
        self.assertAlmostEqual(t.attrs[GTFAttr.ANN_COV_RATIO], 0.2, 2)
        # same strand overlap (no introns)
        t = t_dict['T8']
        self.assertTrue(t.attrs[GTFAttr.CATEGORY] == 0)
        self.assertTrue(t.attrs[GTFAttr.ANN_REF_ID] == 'T1')
        self.assertAlmostEqual(t.attrs[GTFAttr.ANN_INTRON_RATIO], 0.0, 2)
        self.assertAlmostEqual(t.attrs[GTFAttr.ANN_COV_RATIO], 0.10, 2)
        # same strand overlap (with introns)
        t = t_dict['T9']
        self.assertTrue(t.attrs[GTFAttr.CATEGORY] == 0)
        self.assertTrue(t.attrs[GTFAttr.ANN_REF_ID] == 'T1')
        self.assertAlmostEqual(t.attrs[GTFAttr.ANN_INTRON_RATIO], 0.4, 2)
        self.assertAlmostEqual(t.attrs[GTFAttr.ANN_COV_RATIO], 2./3, 2)
        # another same strand overlap (with introns)
        t = t_dict['T10']
        self.assertTrue(t.attrs[GTFAttr.CATEGORY] == 0)
        self.assertTrue(t.attrs[GTFAttr.ANN_REF_ID] == 'T1')
        self.assertAlmostEqual(t.attrs[GTFAttr.ANN_INTRON_RATIO], 0.5, 2)
        self.assertAlmostEqual(t.attrs[GTFAttr.ANN_COV_RATIO], 0.6, 2)

    def test_find_best_match(self):
        transcripts = read_first_locus("annotate_best_match1.gtf")
        t_dict = dict((t.attrs['transcript_id'],t) for t in transcripts)
        annotate_locus(transcripts, 
                       gtf_sample_attr="sample_id") 
        t = t_dict['T1']
        self.assertTrue(t.attrs[GTFAttr.CATEGORY] == 0)
        self.assertTrue(t.attrs[GTFAttr.ANN_REF_ID] == 'D')
        self.assertAlmostEqual(t.attrs[GTFAttr.ANN_INTRON_RATIO], 1.0, 2)
        self.assertAlmostEqual(t.attrs[GTFAttr.ANN_COV_RATIO], 0.9375, 2)
        t = t_dict['T2']
        self.assertTrue(t.attrs[GTFAttr.CATEGORY] == 0)
        self.assertTrue(t.attrs[GTFAttr.ANN_REF_ID] == 'B')
        self.assertAlmostEqual(t.attrs[GTFAttr.ANN_INTRON_RATIO], 0.25, 2)
        self.assertAlmostEqual(t.attrs[GTFAttr.ANN_COV_RATIO], 0.6744, 2)

    def test_intergenic(self):
        transcripts = read_first_locus("annotate_intergenic1.gtf")
        t_dict = dict((t.attrs['transcript_id'],t) for t in transcripts)
        annotate_locus(transcripts, 
                       gtf_sample_attr="sample_id")
        t = t_dict['T1']
        self.assertTrue(t.attrs[GTFAttr.CATEGORY] == 6)
        self.assertTrue(t.attrs[GTFAttr.ANN_REF_ID] == 'na')
        self.assertAlmostEqual(t.attrs[GTFAttr.ANN_INTRON_RATIO], 0.0, 2)
        self.assertAlmostEqual(t.attrs[GTFAttr.ANN_COV_RATIO], 0.0, 2)
        self.assertAlmostEqual(t.attrs[GTFAttr.MEAN_RECURRENCE], 2.0, 2)
        self.assertAlmostEqual(t.attrs[GTFAttr.MEAN_SCORE], 20.0, 2)
        t = t_dict['F']
        self.assertTrue(t.attrs[GTFAttr.CATEGORY] == 6)
        self.assertTrue(t.attrs[GTFAttr.ANN_REF_ID] == 'na')
        self.assertAlmostEqual(t.attrs[GTFAttr.ANN_INTRON_RATIO], 0.0, 2)
        self.assertAlmostEqual(t.attrs[GTFAttr.ANN_COV_RATIO], 0.0, 2)
        self.assertAlmostEqual(t.attrs[GTFAttr.MEAN_RECURRENCE], 2.0, 2)
        self.assertAlmostEqual(t.attrs[GTFAttr.MEAN_SCORE], 20.0, 2)

    def test_test_transcripts(self):
        transcripts = read_first_locus("annotate_test1.gtf")
        t_dict = dict((t.attrs['transcript_id'],t) for t in transcripts)
        annotate_locus(transcripts, 
                       gtf_sample_attr="sample_id")
        t = t_dict['AA']
        self.assertTrue(t.attrs[GTFAttr.CATEGORY] == 1)
        self.assertTrue(t.attrs[GTFAttr.ANN_REF_ID] == "B")
        self.assertTrue(t.attrs[GTFAttr.TEST] == "0")
        t = t_dict['BB']
        self.assertTrue(t.attrs[GTFAttr.CATEGORY] == 0)
        self.assertTrue(t.attrs[GTFAttr.ANN_REF_ID] == "B")
        self.assertTrue(t.attrs[GTFAttr.TEST] == "1")
        t = t_dict['CC']
        self.assertTrue(t.attrs[GTFAttr.CATEGORY] == 1)
        self.assertTrue(t.attrs[GTFAttr.ANN_REF_ID] == "C")
        self.assertTrue(t.attrs[GTFAttr.TEST] == "1")
        t = t_dict['DD']
        self.assertTrue(t.attrs[GTFAttr.CATEGORY] == 3)
        self.assertTrue(t.attrs[GTFAttr.ANN_REF_ID] == "C")
        self.assertTrue(t.attrs[GTFAttr.TEST] == "1")
        t = t_dict['EE']
        self.assertTrue(t.attrs[GTFAttr.CATEGORY] == 5)
        self.assertTrue(t.attrs[GTFAttr.ANN_REF_ID] == "C")
        self.assertTrue(t.attrs[GTFAttr.TEST] == "1")
        t = t_dict['FF']
        self.assertTrue(t.attrs[GTFAttr.CATEGORY] == 6)
        self.assertTrue(t.attrs[GTFAttr.ANN_REF_ID] == "C")
        self.assertTrue(t.attrs[GTFAttr.TEST] == "1")


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()