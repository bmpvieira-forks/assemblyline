'''
Created on Jan 4, 2012

@author: mkiyer
'''
import sys

for line in open(sys.argv[1]):
    fields = line.strip().split('\t')
    start = fields[4]
    if not fields[6]:
        fields[6] = fields[4]
    if not fields[7]:
        fields[7] = fields[4]
    print '\t'.join(fields)