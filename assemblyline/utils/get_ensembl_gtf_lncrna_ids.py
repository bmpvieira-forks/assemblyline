'''
Created on Feb 13, 2013

@author: mkiyer
'''
import argparse

from assemblyline.lib.gtf import GTFFeature

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('gtf_file')
    args = parser.parse_args()
    lncrna_biotypes = set(['3prime_overlapping_ncrna',
                           'antisense', 
                           'lincRNA',
                           'sense_intronic',
                           'sense_overlapping'])
    transcript_ids = set()
    for feature in GTFFeature.parse(open(args.gtf_file)):
        if feature.feature_type != "exon":
            continue
        biotype = feature.attrs["gene_biotype"]
        if biotype in lncrna_biotypes:
            transcript_ids.add(feature.attrs['transcript_id'])
    for t_id in sorted(transcript_ids):
        print t_id
