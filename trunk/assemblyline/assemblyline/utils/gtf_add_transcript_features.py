'''
Created on Jan 22, 2013

@author: mkiyer
'''
import sys
import logging
import operator
import collections
import argparse
import os

from assemblyline.lib.gtf import GTFAttr, GTFFeature, sort_gtf

def gtf_add_transcript_features(gtf_file, outfh):
    transcript_dict = collections.defaultdict(lambda: [])
    for feature in GTFFeature.parse(open(gtf_file)):
        if feature.feature_type != "exon":
            continue
        t_id = feature.attrs[GTFAttr.TRANSCRIPT_ID]
        transcript_dict[t_id].append(feature)
    # output reference transcripts
    for t_id, features in transcript_dict.iteritems():
        # sort features (exons) by start position
        features.sort(key=operator.attrgetter('start'))
        # annotate exons as reference features
        for f in features:
            f.attrs[GTFAttr.REF] = '1'
            print >>outfh, str(f)
        # transcript feature
        f = GTFFeature()
        f.seqid = features[0].seqid
        f.source = features[0].source
        f.feature_type = 'transcript'
        f.start = features[0].start
        f.end = features[-1].end
        f.score = features[0].score
        f.strand = features[0].strand
        f.phase = '.'
        f.attrs = features[0].attrs.copy()
        if "exon_number" in f.attrs:
            del f.attrs["exon_number"]
        f.attrs[GTFAttr.REF] = '1'
        print >>outfh, str(f)

def main():
    # setup logging
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    # parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument('gtf_file')
    args = parser.parse_args()
    # check command line parameters
    if not os.path.exists(args.gtf_file):
        parser.error("GTF file %s not found" % (args.ref_gtf_file))
    # add transcript features
    gtf_add_transcript_features(args.gtf_file, sys.stdout)
    return 0

if __name__ == '__main__':
    sys.exit(main())