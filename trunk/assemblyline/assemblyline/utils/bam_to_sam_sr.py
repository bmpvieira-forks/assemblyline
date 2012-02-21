'''
Created on Feb 20, 2012

@author: mkiyer
'''
import sys
import pysam

if __name__ == '__main__':
    inf = pysam.Samfile(sys.argv[1], "rb")
    outf = pysam.Samfile("-", "w", template=inf)
    for r in inf:
        r.is_paired = False
        r.is_proper_pair = False
        r.is_read1 = True
        outf.write(r)
    inf.close()
    outf.close()
