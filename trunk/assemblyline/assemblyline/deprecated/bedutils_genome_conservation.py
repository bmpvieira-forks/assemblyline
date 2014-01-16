'''
Created on Nov 25, 2013

@author: mkiyer
'''
import logging
import argparse
import os
import sys
import numpy as np

from bedutils_transcript_conservation import which, find_bigwig_files, extract_bigwig_data, BEDFeature, BIGWIG_TO_BEDGRAPH_BIN

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("--pattern", dest="pattern", default=r'{{CHROM}}.phastCons46way.bw')
    parser.add_argument('--bed', dest='bed_file', default=None)
    parser.add_argument("bigwig_file_dir")
    args = parser.parse_args()
    if which(BIGWIG_TO_BEDGRAPH_BIN) is None:
        parser.error('bigWigToBedGraph binary not found in PATH')
    if args.bed_file is not None:
        if not os.path.exists(args.bed_file):
            parser.error('bed file %s not found' % (args.bed_file))
    else:
        parser.error('specify a bed file using --bed')
    # find bigwig files
    chrom_bigwig_dict = find_bigwig_files(args.bigwig_file_dir, args.pattern)
    i = 0
    for f in BEDFeature.parse(open(args.bed_file)):
        # read conservation
        bigwig_file = chrom_bigwig_dict[f.chrom]
        arr = extract_bigwig_data(f, bigwig_file) 
        fields = [f.chrom, str(f.tx_start), str(f.tx_end), str(np.mean(arr))]
        print '\t'.join(fields)
        i += 1
        if (i % 1000) == 0:
            logging.debug("finished %d" % (i))
    return 0

if __name__ == '__main__':
    sys.exit(main())