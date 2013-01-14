'''
Created on Jan 4, 2012

@author: mkiyer
'''
import sys
import logging
import argparse

# setup logging
logging.basicConfig(level=logging.DEBUG,
                    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
# command line parsing
parser = argparse.ArgumentParser()
parser.add_argument("chrom_sizes_file")
parser.add_argument("rmsk_file")
args = parser.parse_args()
chrominfo = args.chrom_sizes_file
rmsk = args.rmsk_file

chromsizes = {}
for line in open(chrominfo):
    fields = line.strip().split()
    chrom = fields[0]
    length = int(fields[1])
    chromsizes[chrom] = length
    print "@SQ\tSN:%s\tLN:%d" % (chrom, length)

for line in open(rmsk):
    if line.startswith("#"):
        continue
    fields = line.strip().split()
    chrom = fields[5]
    start = int(fields[6])
    end = int(fields[7])
    strand = fields[9]
    repname = fields[10]
    repclass = fields[11]
    repfamily = fields[12]
    if repfamily != "rRNA":
        continue
    name = "%s:%d-%d_%s_%s" % (chrom, start, end, strand, repname)
    if chrom not in chromsizes:
        print >>sys.stderr, "could not find chrom", chrom
        continue
    print '\t'.join(map(str, [chrom, start, end, strand, name]))