'''
Created on Nov 25, 2013

@author: mkiyer
'''
import logging
import argparse
import os
import sys
import random
import bisect

from base import BEDFeature

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', type=int, dest='num_samples', default=100000)
    parser.add_argument('genome_bed_file')
    parser.add_argument('region_lengths_file')
    # parse arguments
    args = parser.parse_args()
    genome_bed_file = args.genome_bed_file
    region_lengths_file = args.region_lengths_file
    num_samples = max(1, args.num_samples)
    # check arguments
    if not os.path.exists(genome_bed_file):
        parser.error('genome bed file %s not found' % (genome_bed_file))
    if not os.path.exists(region_lengths_file):
        parser.error('region lengths file %s not found' % (region_lengths_file))
    # get transcript lengths from transcript bed file
    logging.debug('Reading region lengths')
    region_lengths = []
    with open(region_lengths_file) as f:
        region_lengths.extend(int(x.strip()) for x in f)
    # first pass to index windows
    logging.debug('Reading genome bed file')
    features = []
    genome_cumsums = []
    genome_size = 0
    max_feature_size = 0
    for f in BEDFeature.parse(open(genome_bed_file)):
        length = sum((e[1]-e[0]) for e in f.exons)
        if length == 0:
            continue
        features.append(f)
        genome_cumsums.append(genome_size)
        genome_size += length
        max_feature_size = max(max_feature_size, length)
    logging.debug('Genome bed size %d' % (genome_size))
    # get windows
    trials = 0
    windows = []
    for i in xrange(num_samples):
        # choose random region length
        length = region_lengths[random.randrange(len(region_lengths))]        
        while True:
            trials += 1
            # choose random window start                
            genome_start = random.randrange(genome_size - length)        
            feature_index = bisect.bisect_right(genome_cumsums, genome_start) - 1
            f = features[feature_index]
            feature_length = (f.tx_end - f.tx_start)
            # can use available length of this feature
            genome_feature_start = genome_cumsums[feature_index]
            feature_offset = genome_start - genome_feature_start
            avail_length = feature_length - feature_offset
            if avail_length < length:
                continue
            windows.append((f.chrom, f.tx_start + feature_offset, f.tx_start + feature_offset + length))
            break
        if (i > 0) and (i % 1000) == 0:
            logging.debug('Finished %d' % (i))
    logging.debug('Sampled %d windows (%d trials)' % (len(windows), trials))
    for window in sorted(windows):
        print '\t'.join(map(str, window))
    return 0

if __name__ == '__main__':
    sys.exit(main())