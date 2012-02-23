'''
Created on Feb 20, 2012

@author: mkiyer
'''
import sys
import pysam

if __name__ == '__main__':
    fin = pysam.Samfile(sys.argv[1], "rb")
    fout = pysam.Samfile("-", "w", template=fin)
    for r in fin:
        if r.is_paired:
            r.is_paired = False
            r.is_proper_pair = False
            r.is_read1 = True
            fout.write(r)
    fin.close()
    fout.close()
