'''
Created on Feb 15, 2011

@author: mkiyer
'''
import os

def up_to_date(outfile, infile, nzsize=True):
    if not os.path.exists(infile):
        return False
    if not os.path.exists(outfile):
        return False
    if nzsize and (os.path.getsize(outfile) == 0):
        return False
    return os.path.getmtime(outfile) >= os.path.getmtime(infile)

imax2 = lambda x,y: x if x>=y else y
imin2 = lambda x,y: x if x<=y else y

class GTFAttr:
    GENE_ID = 'gene_id'
    TRANSCRIPT_ID = 'transcript_id'
    COHORT_ID = 'cohort_id'
    SAMPLE_ID = 'sample_id'
    LIBRARY_ID = 'library_id'
    REF = 'ref'

