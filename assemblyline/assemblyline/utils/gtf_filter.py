'''
Created on Dec 31, 2013

@author: mkiyer
'''
import argparse
import sys

chrom_names = set(["chr1", "1",
                   "chr2", "2",
                   "chr3", "3",
                   "chr4", "4",
                   "chr5", "5",
                   "chr6", "6",
                   "chr7", "7",
                   "chr8", "8",
                   "chr9", "9",
                   "chr10", "10",
                   "chr11", "11",
                   "chr12", "12",
                   "chr13", "13",
                   "chr14", "14",
                   "chr15", "15",
                   "chr16", "16",
                   "chr17", "17",
                   "chr18", "18",
                   "chr19", "19",
                   "chr20", "20",
                   "chr21", "21",
                   "chr22", "22",
                   "chrX", "X",
                   "chrY", "Y",
                   "chrM", "MT"])

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("filename")
    args = parser.parse_args()
    skipped = 0
    with open(args.filename) as f:
        for line in f:
            fields = line.strip().split('\t')
            if fields[0] in chrom_names:
                print '\t'.join(fields)
            else:
                skipped += 1
    print >>sys.stderr, 'skipped %d lines' % (skipped)
    
if __name__ == '__main__':
    sys.exit(main())
