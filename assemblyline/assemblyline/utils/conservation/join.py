'''
Created on Jan 8, 2014

@author: mkiyer
'''
import sys
import argparse

parser = argparse.ArgumentParser()

grp1 = parser.add_mutually_exclusive_group(required=True)
grp1.add_argument('-1', '--f1', dest='field1', default=None)
grp1.add_argument('-n1', dest='col1', type=int, default=None)
grp2 = parser.add_mutually_exclusive_group(required=True)
grp2.add_argument('-2', '-f2', dest='field2', default=None)
grp2.add_argument('-n2', dest='col2', type=int, default=None)
parser.add_argument('--sep', dest='sep', default='\t')
parser.add_argument('file1')
parser.add_argument('file2')
args = parser.parse_args()

if args.col1 is not None:
    index1 = args.col1
else:
    index1 = None
if args.col2 is not None:
    index2 = args.col2
else:
    index2 = None
file1 = args.file1
file2 = args.file2
field1 = args.field1
field2 = args.field2
sep = args.sep

d = {}
header2 = []
numfields2 = None
with open(file2) as f:
    if index2 is None:
        header2 = f.next().strip().split(sep)
        index2 = header2.index(field2)
        del header2[index2]
        numfields2 = len(header2)
    for line in f:
        fields = line.strip().split(sep)
        key = fields[index2]
        del fields[index2]
        d[key] = fields
        if numfields2 is None:
            numfields2 = len(fields)

header1 = []
missing = 0
blank = ['NA'] * numfields2
with open(file1) as f:
    if index1 is None:
        header1 = f.next().strip().split(sep)
        index1 = header1.index(field1)
    if len(header1) > 0 and len(header2) > 0:
        print sep.join(header1 + header2)
    for line in f:
        fields = line.strip().split(sep)
        key = fields[index1]
        if key not in d:
            fields.extend(blank)
            missing += 1
        else:
            fields.extend(d[key])
        print sep.join(fields)

print >>sys.stderr, "found %d missing keys" % (missing)