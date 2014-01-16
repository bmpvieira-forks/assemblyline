'''
Created on Jan 7, 2014

@author: mkiyer
'''
import argparse
import logging
import collections

import numpy as np
import matplotlib.pyplot as plt

from conservation_histogram import NUM_BINS, BINS, BIN_MIN, BIN_MAX

def get_hists_promoters(npzfile, hists):
    for key in npzfile.keys():
        hist = np.array(npzfile[key])
        groups = key.split('|')
        if groups[0] == 'random':
            hists[groups[2]] += hist
        elif len(groups) == 1:
            tcat, cat = groups[0].split('.')
            if tcat == 'lncRNA':
                if cat == 'same_strand':
                    hists['known lncRNA'] += hist
                else:
                    hists['novel lncRNA'] += hist
            elif tcat == 'protein_coding':
                hists[tcat] += hist

def get_hists_transcripts(npzfile, hists):
    for key in npzfile.keys():
        print key
        continue
        tcat, mcat, ccat = key.split('|')
        hist = np.array(npzfile[key])
        #hists[key] += hist
        #continue
        if (ccat == 'noncoding'):
            if tcat == 'lncRNA':
                k = tcat
                #if mcat == 'same_strand':
                #    k = 'known lncRNA'
                #else:
                #    k = 'novel lncRNA'
            elif tcat == 'pseudogene':
                k = tcat
            else:
                k = 'other ncRNA'
        else:
            k = ccat
        hists[k] += hist


def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    # parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("histogram_files", nargs="+")
    args = parser.parse_args()
    filenames = args.histogram_files
    hists = collections.defaultdict(lambda: np.zeros(NUM_BINS-1, dtype=np.float))
    for filename in filenames:
        npzfile = np.load(filename)
        get_hists_promoters(npzfile, hists)
    # plot
    for k in sorted(hists):
        hist = hists[k]
        density = hist.cumsum() / hist.sum()
        plt.plot(BINS[:-1], density, label=k)
        print k, hist.sum()
    plt.xlim((-4.0,6.0))
    plt.legend(loc="upper left")
    plt.show()
        

if __name__ == '__main__':
    main()
