'''
Created on Feb 20, 2012

@author: mkiyer
'''
import sys
import pysam

if __name__ == '__main__':
    fin = pysam.Samfile(sys.argv[1], "rb")
    fout = pysam.Samfile(sys.argv[2], "wh", template=fin)
    for r in fin:
        suffix = r.qname[-2:]
        if (suffix == "/1" or suffix == "/2"):
            r.qname = r.qname[:-2]
        fout.write(r)
    fin.close()
    fout.close()
