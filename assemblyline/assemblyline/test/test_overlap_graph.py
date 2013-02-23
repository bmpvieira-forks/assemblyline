'''
Created on Feb 22, 2013

@author: mkiyer
'''
import unittest

from assemblyline.lib.transcript import Exon
from assemblyline.lib.assemble.assembler import get_nmers

class TestOverlapGraph(unittest.TestCase):

    def test_get_nmers(self):
        
        path = (Exon(0,10),Exon(10,15),Exon(15,20),Exon(20,70),
        
        pass


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()