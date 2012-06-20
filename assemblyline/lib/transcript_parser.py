'''
Created on Nov 2, 2010

@author: mkiyer
'''
import logging
import operator
import collections
import os
import math

import gtf
from transcript import Exon, Transcript, strand_str_to_int

def float_check_nan(x):
    x = float(x)
    if math.isnan(x):
        return 0.0
    return x

cufflinks_attr_defs = {"FPKM": float_check_nan,
                       "cov": float_check_nan,
                       "frac": float_check_nan}

def separate_transcripts(gtf_features):
    transcripts = collections.defaultdict(lambda: Transcript())
    for feature in gtf_features:
        tx_id = feature.attrs["transcript_id"]
        if feature.feature_type == "exon" and tx_id not in transcripts:            
            logging.error("Feature type 'exon' has no transcript record: %s" % feature)
            assert False
        transcript = transcripts[tx_id]
        if feature.feature_type == "transcript":
            transcript.chrom = feature.seqid
            transcript.start = feature.start
            transcript.end = feature.end
            # convert from string strand notation ("+", "-", ".") 
            # to integer (0, 1)            
            transcript.strand = strand_str_to_int(feature.strand)
            transcript.exons = {}
            transcript.attrs = feature.attrs
        elif feature.feature_type == "exon":
            exon_num = int(feature.attrs["exon_number"])
            transcript.exons[exon_num] = Exon(feature.start, feature.end)
    # convert transcript exons from dictionary to sorted list
    for transcript in transcripts.itervalues():
        # sort the exons in correct transcript order
        transcript.exons = [transcript.exons[v] for v in sorted(transcript.exons.keys())]
    # sort transcripts by position
    return sorted(transcripts.values(), key=operator.attrgetter('start'))

def window_overlap(a, b):
    if a[0] != b[0]:
        return False
    return (a[1] <= b[2]) and (b[1] <= a[2])

def separate_loci(feature_iter):
    try:
        # initialize window
        window = [feature_iter.next()]
        window_range = (window[0].seqid, window[0].start, window[0].end)
        # separate into loci
        for feature in feature_iter:
            # check if next transcript is outside current window
            interval = (feature.seqid, feature.start, feature.end)
            if not window_overlap(interval, window_range):
                # yield current window
                yield window
                # reset window
                window = [feature]
                window_range = (feature.seqid, feature.start, feature.end)
            else:
                # add transcript to window
                window.append(feature)
                window_range = (feature.seqid,
                                min(window_range[1], feature.start),
                                max(window_range[2], feature.end))
    except StopIteration:
        pass
    # yield last window
    if len(window) > 0:
        yield window

def parse_gtf(fileh, attr_defs=None):
    for locus_features in separate_loci(gtf.GTFFeature.parse(fileh, attr_defs)):
        locus_transcripts = separate_transcripts(locus_features)
        yield locus_transcripts

if __name__ == '__main__':
    import argparse
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")    
    parser = argparse.ArgumentParser()
    parser.add_argument("filename")
    options = parser.parse_args()
    prefix = os.path.splitext(options.filename)[0]
    output_file = options.filename
    for locus_transcripts in parse_gtf(open(output_file), attr_defs=None):        
        for t in locus_transcripts:
            print t
        print '------------'
