'''
Created on Feb 23, 2013

@author: mkiyer
'''
import argparse

from assemblyline.lib.gtf import GTFFeature

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('gtf_file')
    args = parser.parse_args()
    transcript_ids = set()
    for feature in GTFFeature.parse(open(args.gtf_file)):
        if feature.feature_type != "exon":
            continue
        transcript_ids.add(feature.attrs['transcript_id'])
    for t_id in sorted(transcript_ids):
        print t_id
