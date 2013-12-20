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
import numpy as np

from bedutils_transcript_conservation import extract_bigwig_data, BEDFeature

def find_valid_windows(arr, window_size, threshold, score_func=np.mean):
    windows = []
    if window_size > len(arr):
        return windows
    if len(arr) == 0:
        return windows
    num_windows = len(arr) - window_size + 1
    for i in xrange(num_windows):
        score = score_func(arr[i:i+window_size])
        if score >= threshold:
            windows.append(i)
    return windows

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("--bin-path", dest="bin_path", default='')
    parser.add_argument('--window-size', type=int, dest='window_size', default=300)
    parser.add_argument('--threshold', type=float, dest='mapability_threshold', default=0.99)
    parser.add_argument('-n', type=int, dest='num_samples', default=1000)
    parser.add_argument('bigwig_file')
    parser.add_argument("bed_file")
    args = parser.parse_args()
    bigwig_to_bedgraph_bin = os.path.join(args.bin_path, 'bigWigToBedGraph')
    if not os.path.exists(bigwig_to_bedgraph_bin):
        parser.error('bigwig to bedgraph converter binary not found')
    if not os.path.exists(args.bed_file):
        parser.error('bed file %s not found' % (args.bed_file))
    if not os.path.exists(args.bigwig_file):
        parser.error('bigwig file %s not found' % (args.bigwig_file))
    bigwig_file = args.bigwig_file
    window_size = max(1, args.window_size)
    mapability_threshold = args.mapability_threshold
    num_samples = max(1, args.num_samples)
    # parameters
    logging.info("Window size: %d" % (window_size))
    logging.info("Mapability file: %s" % (bigwig_file))
    logging.info("Mapability threshold: %f" % (mapability_threshold))
    # first pass to index windows
    features = []
    window_cumsums = []
    num_windows = 0
    for f in BEDFeature.parse(open(args.bed_file)):
        length = sum((e[1]-e[0]) for e in f.exons)
        if length == 0:
            continue
        if window_size > length:
            continue
        window_cumsums.append(num_windows)
        features.append(f)
        num_windows += length - window_size + 1
    logging.debug('Found %d total windows' % (num_windows))
    # now sample windows
    if num_windows < num_samples:
        logging.warning('Only sampling %d windows' % (num_windows))
    num_samples = min(num_samples, num_windows)
    windows = []
    trials = 0
    while len(windows) < num_samples:
        # choose random window
        windex = random.randrange(num_windows)
        feature_index = bisect.bisect_right(window_cumsums, windex) - 1
        woffset = windex - window_cumsums[feature_index]
        f = features[feature_index]
        # fetch data in window
        arr = extract_bigwig_data(f, bigwig_file, bigwig_to_bedgraph_bin) 
        score = np.mean(arr[woffset:woffset+window_size])
        if score >= mapability_threshold:
            windows.append((f.chrom, f.tx_start + woffset, f.tx_start + woffset + window_size))
            logging.debug('Sampled %d windows (%d trials)' % (len(windows), trials))
        trials += 1
    for window in sorted(windows):
        print '\t'.join(map(str, window))
    return 0

if __name__ == '__main__':
    sys.exit(main())