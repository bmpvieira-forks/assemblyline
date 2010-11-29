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

def merge_strand(strand1, strand2):
    if strand1 == strand2:
        return strand1
    elif strand1 == ".":
        return strand2
    elif strand2 == ".":
        return strand1
    else:
        logging.error("Incompatible strands")        
        assert False

class Transcript(object):
    __slots__ = ('chrom', 'start', 'end', 'strand', 'id', 'label', 'score', 'exons')
    def __init__(self):
        self.id = None
        self.label = None
        self.chrom = None
        self.start = -1
        self.end = -1
        self.strand = '.'
        self.score = 0
        self.exons = {}
    def __str__(self):
        return ("<%s(chrom='%s', start=%d, end=%d, strand=%s, id=%s, label=%s, score=%f, exons=%s)>" %
                (self.__class__.__name__, self.chrom, self.start, self.end, self.strand, self.id, self.label, self.score, self.exons))
    @property
    def tx_start(self):
        return self.exons[0][0]
    @property
    def tx_end(self):
        return self.exons[-1][1]
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
                       self.strand if self.strand != None else '.',
                       str(self.tx_start),
                       str(self.tx_start),
                       '0',
                       str(len(self.exons)),
                       ','.join(map(str,block_sizes)) + ',',
                       ','.join(map(str,block_starts)) + ','])
        return s



class Exon(object):
    __slots__ = ('start', 'end', 'score')
    def __str__(self):
        return ("<%s(start=%d, end=%d, score=%f)>" %
                (self.__class__.__name__, self.start, self.end, self.score))

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
            transcript.strand = feature.strand
            transcript.chrom = feature.seqid
            transcript.start = feature.start
            transcript.end = feature.end
            transcript.score = float(feature.attrs[score_attr])
        elif feature.feature_type == "exon":
            exon_num = int(feature.attrs["exon_number"])
            exon = Exon()
            exon.start = feature.start
            exon.end = feature.end
            # remove the length normalization from scores
            exon.score = float(feature.attrs[score_attr]) * (feature.end - feature.start) / 1.0e3
            transcript.exons[exon_num] = exon
    # convert transcript exons from dictionary to sorted list
    for transcript in transcripts.itervalues():
        # sort the exons in correct transcript order
        transcript.exons = [transcript.exons[v] for v in sorted(transcript.exons.keys())]
        # remove the length normalization from scores
        transcript_length = sum([e.end - e.start for e in transcript.exons])
        transcript.score = transcript.score * (transcript_length) / 1.0e3
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
    logging.debug("Reading features")
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
