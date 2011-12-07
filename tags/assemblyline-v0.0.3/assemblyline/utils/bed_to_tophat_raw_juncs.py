'''
Created on Oct 26, 2011

@author: mkiyer
'''
import os
import sys
import argparse
import logging

def bed_to_tophat_juncs(line, linenum):
    if not line:
        return
    if line.startswith("#") or line.startswith("track"):
        return
    fields = line.strip().split('\t')
    if len(fields) < 12:
        logging.warning("malformed line %d, missing columns" % linenum)
    chrom = fields[0]
    start = int(fields[1])
    strand = fields[5]
    block_sizes = map(int, fields[10].split(','))
    block_starts = map(int, fields[11].split(','))
    for i in xrange(1,len(block_sizes)):
        if not block_sizes[i]:
            continue
        left_pos = start + block_starts[i-1] + block_sizes[i-1] - 1
        right_pos = start + block_starts[i]
        yield chrom, left_pos, right_pos, strand

def get_refs_from_bowtie_index(bowtie_index, bowtie_inspect_bin='bowtie-inspect'):
    import subprocess
    args = [bowtie_inspect_bin, '-s', bowtie_index]    
    p = subprocess.Popen(args, stdout=subprocess.PIPE)
    output = p.communicate()[0]
    refs = []
    for line in output.split('\n'):
        if not line:
            continue
        fields = line.split('\t')
        if fields[0].startswith('Sequence'):
            refs.append((fields[1], int(fields[2])))
    return refs

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', dest="output_file", default=None)
    parser.add_argument('--bowtie-inspect-bin', dest="bowtie_inspect_bin", default="bowtie-inspect")
    parser.add_argument('index')    
    parser.add_argument('bed_file')
    args = parser.parse_args()     
    # check args
    if not os.path.exists(args.index + ".1.ebwt"):
        parser.error("bowtie index %s not found.. exiting" % (args.index))    
    if not os.path.exists(args.bed_file):
        parser.error("BED file %s not found... exiting." % (args.bed_file))
    # getting references
    refs = get_refs_from_bowtie_index(args.index, bowtie_inspect_bin=args.bowtie_inspect_bin)
    refs = set(x[0] for x in refs)
    if args.output_file is not None:
        f = open(args.output_file, "w")
    else:
        f = sys.stdout
    juncs = set()
    for linenum,line in enumerate(open(args.bed_file)):
        for chrom, left, right, strand in bed_to_tophat_juncs(line, linenum):
            if chrom not in refs:
                logging.warning("skipping line %d: %s" % (linenum, line.strip()))
            juncs.add((chrom, left, right, strand))
    for chrom, left, right, strand in sorted(juncs):
        print >>f, "%s\t%d\t%d\t%s" % (chrom, left, right, strand)
    if args.output_file is not None:
        f.close()

if __name__ == '__main__':
    main()
    #sys.exit(main())