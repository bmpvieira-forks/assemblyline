'''
Created on Jan 30, 2014

@author: mkiyer
'''
import logging
import argparse
import os
import sys

from assemblyline.utils.conservation.base import BEDFeature
from assemblyline.lib.gtf import GTFFeature

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", dest="source", default='bed_to_gtf')
    parser.add_argument("bed_file")
    args = parser.parse_args()
    bed_file = args.bed_file
    source = args.source
    f_num = 1
    for x in BEDFeature.parse(open(bed_file)):
        f = GTFFeature()
        f.seqid = x.chrom
        f.source = source
        f.feature_type = 'transcript'
        f.start = x.tx_start
        f.end = x.tx_end
        f.score = x.score
        f.strand = x.strand
        f.phase = '.'
        f.attrs = {'transcript_id': 'T%08d' % (f_num),
                   'gene_id': 'G%08d' % (f_num)}
        features = [f]
        for i,e in enumerate(x.exons):
            start, end = e
            f = GTFFeature()
            f.seqid = x.chrom
            f.source = source
            f.feature_type = 'exon'
            f.start = start
            f.end = end
            f.score = x.score
            f.strand = x.strand
            f.phase = '.'
            f.attrs = dict(features[0].attrs)
            f.attrs["exon_number"] = i
            features.append(f)
        for f in features:
            print str(f)
        f_num += 1

if __name__ == '__main__':
    sys.exit(main())
            
