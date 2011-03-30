'''
Created on Nov 2, 2010

@author: mkiyer
'''
import logging
import operator
import collections
import subprocess
import os
import sys

import gtf
from base import Exon, NO_STRAND, POS_STRAND, NEG_STRAND, strand_int_to_str, \
    strand_str_to_int

class Transcript(object):
    __slots__ = ('chrom', 'start', 'end', 'strand', 'id', 'label', 
                 'score', 'length', 'exons', 'attrs')
    def __init__(self):
        self.chrom = None
        self.start = -1
        self.end = -1
        self.strand = NO_STRAND
        self.id = None
        self.label = None
        self.score = 0
        self.length = 0
        self.exons = {}
        self.attrs = {}
    def __str__(self):
        return ("<%s(chrom='%s', start='%d', end='%d', strand='%s', "
                "id='%s', label='%s', score='%f', length='%d', exons=%s, attrs='%s')>" %
                (self.__class__.__name__, self.chrom, self.start, self.end, 
                 strand_int_to_str(self.strand), self.id, self.label, 
                 self.score, self.length, self.exons, self.attrs))
    @property
    def tx_start(self):
        return self.exons[0].start
    @property
    def tx_end(self):
        return self.exons[-1].end
    @property
    def introns(self):
        return zip([e[1] for e in self.exons[:-1]],
                   [e[0] for e in self.exons[1:]])

    def to_bed12(self):
        block_sizes = []
        block_starts = []
        for e0, e1 in self.exons:
            block_starts.append(e0 - self.tx_start)
            block_sizes.append(e1 - e0)        
        # write
        s = '\t'.join([self.chrom, 
                       str(self.tx_start), 
                       str(self.tx_end),
                       str(self.id),
                       str(self.score),
                       strand_int_to_str(self.strand),
                       str(self.tx_start),
                       str(self.tx_start),
                       '0',
                       str(len(self.exons)),
                       ','.join(map(str,block_sizes)) + ',',
                       ','.join(map(str,block_starts)) + ','])
        return s


def separate_transcripts(gtf_features, score_attr="FPKM"):
    transcripts = collections.defaultdict(lambda: Transcript())
    for feature in gtf_features:
        tx_id = feature.attrs["transcript_id"]
        seqdata_id = feature.attrs["gene_id"].split(".")[0]
        if feature.feature_type == "exon" and tx_id not in transcripts:            
            logging.error("Feature type 'exon' has no transcript record: %s" % feature)
            assert False
        transcript = transcripts[tx_id]
        if feature.feature_type == "transcript":
            transcript.id = tx_id
            transcript.label = seqdata_id
            # convert from string strand notation ("+", "-", ".") 
            # to integer (0, 1)            
            transcript.strand = strand_str_to_int(feature.strand)
            transcript.chrom = feature.seqid
            transcript.start = feature.start
            transcript.end = feature.end
            transcript.score = float(feature.attrs[score_attr])
        elif feature.feature_type == "exon":
            exon_num = int(feature.attrs["exon_number"])
            transcript.exons[exon_num] = Exon(feature.start, feature.end)
    # convert transcript exons from dictionary to sorted list
    for transcript in transcripts.itervalues():
        # sort the exons in correct transcript order
        transcript.exons = [transcript.exons[v] for v in sorted(transcript.exons.keys())]
        # remove the length normalization from scores
        transcript.length = sum([e.end - e.start for e in transcript.exons])
        transcript.score = transcript.score * (transcript.length) / 1.0e3        
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

def parse_gtf(fileh, attr_defs=None, score_attr="FPKM"):
    for locus_features in separate_loci(gtf.GTFFeature.parse(fileh, attr_defs)):
        locus_transcripts = separate_transcripts(locus_features, score_attr)
        yield locus_transcripts

def sort_gtf(filename, output_file):
    args = ["sort", "-k1,1", "-k4,4n", "-k3,3r", filename]
    myenv = os.environ.copy()
    myenv["LC_ALL"] = "C"
    subprocess.call(args, stdout=open(output_file, "w"), env=myenv)

if __name__ == '__main__':
    import argparse
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")    
    parser = argparse.ArgumentParser()
    parser.add_argument("--sort", action="store_true", default=False)
    parser.add_argument("--score-attr", dest="score_attr", default="FPKM")
    parser.add_argument("filename")
    options = parser.parse_args()
    prefix = os.path.splitext(options.filename)[0]
    if options.sort:
        output_file = prefix + ".srt.gtf"
        sort_gtf(options.filename, output_file)
    else:
        output_file = options.filename
    for locus_transcripts in parse_gtf(open(output_file), attr_defs=None, score_attr=options.score_attr):        
        for t in locus_transcripts:
            print t
        print '------------'
