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

from bedutils_transcript_conservation import extract_bigwig_data, BEDFeature, which, BIGWIG_TO_BEDGRAPH_BIN

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
    parser.add_argument('genome_bed_file')
    parser.add_argument('transcript_bed_file')
    # parse arguments
    args = parser.parse_args()
    bigwig_file = args.bigwig_file
    genome_bed_file = args.genome_bed_file
    transcript_bed_file = args.transcript_bed_file
    window_size = max(1, args.window_size)
    mapability_threshold = args.mapability_threshold
    num_samples = max(1, args.num_samples)
    # check arguments
    if which(BIGWIG_TO_BEDGRAPH_BIN) is None:
        parser.error('bigWigToBedGraph binary not found in PATH')
    if not os.path.exists(bigwig_file):
        parser.error('bigwig file %s not found' % (bigwig_file))
    if not os.path.exists(genome_bed_file):
        parser.error('genome bed file %s not found' % (genome_bed_file))
    if not os.path.exists(transcript_bed_file):
        parser.error('transcript bed file %s not found' % (transcript_bed_file))
    # parameters
    logging.info("Window size: %d" % (window_size))
    logging.info("Mapability file: %s" % (bigwig_file))
    logging.info("Mapability threshold: %f" % (mapability_threshold))
    # get transcript lengths from transcript bed file
    transcript_lengths = []
    for f in BEDFeature.parse(open(transcript_bed_file)):
        length = sum((e[1]-e[0]) for e in f.exons)
        transcript_lengths.append(length)
    # first pass to index windows
    features = []
    window_cumsums = []
    num_windows = 0
    for f in BEDFeature.parse(open(genome_bed_file)):
        length = sum((e[1]-e[0]) for e in f.exons)
        if length == 0:
            continue
        if window_size > length:
            continue
        window_cumsums.append(num_windows)
        features.append(f)
        num_windows += length - window_size + 1
    logging.debug('Found %d total windows' % (num_windows))

    # get windows
    for i in xrange(num_samples):
        
    
    
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
        arr = extract_bigwig_data(f, bigwig_file) 
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