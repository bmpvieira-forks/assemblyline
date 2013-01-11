'''
Created on Jan 10, 2013

@author: mkiyer
'''
import sys
import os
import logging
import argparse
import array

# local imports
import pysam

# project imports
from assemblyline.rnaseq.lib.base import parse_sam, \
    remove_multihits, to_fastq

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("sam_file")
    parser.add_argument("bam_file")
    parser.add_argument("counts_file")
    parser.add_argument("fastq_files")
    args = parser.parse_args()
    # open input files
    suffix = os.path.splitext(args.sam_file)[-1]
    mode = 'rb' if suffix == '.bam' else 'r'
    infh = pysam.Samfile(args.sam_file, mode)
    nrefs = len(infh.references)
    # open output files
    bamfh = pysam.Samfile(args.bam_file, "wb", template=infh)
    fastq_files = args.fastq_files.split(',')
    fastq_fhs = [open(f,'w') for f in fastq_files]
    # init counts
    def make_array(n):
        return array.array('L', (0 for x in xrange(n)))
    counts = {'UU': make_array(nrefs),
              'CP': make_array(nrefs),
              'DP': make_array(nrefs),
              'UP': make_array(nrefs)}
    both_unmapped = 0
    for pe_reads in parse_sam(infh, 
                              readnum_in_qname=False, 
                              remove_suffix=True):
        r1, r2 = remove_multihits(pe_reads)
        assert (r1 is not None) and (r2 is not None)
        # if both reads are unmapped then write to fastq
        if (r1.is_unmapped and r2.is_unmapped):
            print >>fastq_fhs[0], to_fastq(r1, 0)
            print >>fastq_fhs[1], to_fastq(r2, 1)
            both_unmapped += 2
        else:
            for r in (r1, r2):
                if not r.is_unmapped:
                    # 'YT' is bowtie2 SAM tag that gives PE alignment 
                    # information
                    altype = r.opt('YT')
                    counts[altype][r.tid] += 1            
                bamfh.write(r)
    # cleanup
    for f in fastq_fhs:
        f.close()
    bamfh.close()
    # print output
    outfh = open(args.counts_file, "w")
    print >>outfh, '\t'.join(["#rname", "length", "paired", 
                              "paired_discordant", "paired_orphan", 
                              "unpaired", "total"])
    for i in xrange(len(infh.references)):
        paired = counts['CP'][i]
        paired_discordant = counts['DP'][i]
        paired_orphan = counts['UP'][i]
        unpaired = counts['UU'][i]
        total = paired + paired_discordant + paired_orphan + unpaired
        fields = [infh.references[i],
                  infh.lengths[i],
                  paired,
                  paired_discordant,
                  paired_orphan,
                  unpaired,
                  total]
        print >>outfh, '\t'.join(map(str, fields))
    print >>outfh, '\t'.join(['*', '0', '0', '0', '0', '0', 
                              str(both_unmapped)])
    outfh.close()
    infh.close()
    return 0

if __name__ == '__main__':
    sys.exit(main())