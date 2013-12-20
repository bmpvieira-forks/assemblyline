'''
Created on Nov 25, 2013

@author: mkiyer
'''
import logging
import argparse
import os
import sys
import numpy as np

from bedutils_transcript_conservation import find_bigwig_files, extract_bigwig_data, BEDFeature

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("--bin-path", dest="bin_path", default='')
    parser.add_argument("--pattern", dest="pattern", default=r'${CHROM}.phastCons46way.bw')
    parser.add_argument("bigwig_file_dir")
    parser.add_argument("bed_file")
    args = parser.parse_args()
    bigwig_to_bedgraph_bin = os.path.join(args.bin_path, 'bigWigToBedGraph')
    if not os.path.exists(bigwig_to_bedgraph_bin):
        parser.error('bigwig to bedgraph converter binary not found')
    if not os.path.exists(args.bed_file):
        parser.error('bed file %s not found' % (args.bed_file))
    # find bigwig files
    chrom_bigwig_dict = find_bigwig_files(args.bigwig_file_dir, args.pattern)
    for f in BEDFeature.parse(open(args.bed_file)):
        # read conservation
        bigwig_file = chrom_bigwig_dict[f.chrom]
        arr = extract_bigwig_data(f, bigwig_file, bigwig_to_bedgraph_bin) 
        fields = [f.chrom, str(f.tx_start), str(f.tx_end), str(np.mean(arr))]
        print '\t'.join(fields)
    return 0

if __name__ == '__main__':
    sys.exit(main())