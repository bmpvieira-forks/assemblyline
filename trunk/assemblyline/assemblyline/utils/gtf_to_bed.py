'''
Created on Jul 31, 2013

@author: mkiyer
'''
import logging
import argparse
import os
import sys

# project imports
from assemblyline.lib.gtf import GTFFeature
from assemblyline.lib.transcript import parse_gtf, strand_int_to_str

def write_bed(chrom, name, strand, score, exons):
    assert all(exons[0].start < x.start for x in exons[1:])
    assert all(exons[-1].end > x.end for x in exons[:-1])
    tx_start = exons[0].start
    tx_end = exons[-1].end    
    block_sizes = []
    block_starts = []
    for e in exons:
        block_starts.append(e.start - tx_start)
        block_sizes.append(e.end - e.start)        
    # make bed fields
    fields = [chrom, 
              str(tx_start), 
              str(tx_end),
              str(name),
              str(score),
              strand_int_to_str(strand),
              str(tx_start),
              str(tx_start),
              '0',
              str(len(exons)),
              ','.join(map(str,block_sizes)) + ',',
              ','.join(map(str,block_starts)) + ',']
    return fields

def main():
    # setup logging
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    # parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument('-g', dest='gene_id', action='store_true', default=False)
    parser.add_argument('gtf_file')
    args = parser.parse_args()
    # check command line parameters
    if not os.path.exists(args.gtf_file):
        parser.error("gtf file %s not found" % (args.gtf_file))
    for transcripts in parse_gtf(open(args.gtf_file)):
        for t in transcripts:
            if args.gene_id:
                name = '%s|%s' % (t.attrs['gene_id'], t.attrs['transcript_id'])
            else:
                name = t.attrs['transcript_id']
            fields = write_bed(t.chrom, name, t.strand, 1000, t.exons)
            print '\t'.join(fields)
    return 0

if __name__ == '__main__':
    sys.exit(main())
