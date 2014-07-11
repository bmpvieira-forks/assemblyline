'''
Created on Feb 24, 2013

@author: mkiyer
'''
import logging
import argparse
import os
import sys

from assemblyline.lib.gtf import GTFFeature

FULL_GTF_ATTRS = ['gene_id', 
                  'tss_id', 
                  'transcript_name', 
                  'transcript_category', 
                  'transcript_type', 
                  'ref_transcript_id', 
                  'ref_orig_gene_id', 
                  'ref_locus',
                  'ref_source',
                  'ref_gene_name',
                  'ref_gene_type',
                  'ref_length',
                  'ref_num_introns',
                  'category',
                  'distance',
                  'shared_same_strand_bp',
                  'shared_opp_strand_bp',
                  'shared_introns',
                  'shared_splicing']

DEFAULT_GTF_ATTRS = ['gene_id', 
                     'tss_id', 
                     'transcript_name', 
                     'transcript_category', 
                     'transcript_type', 
                     'category',
                     'distance']

class TranscriptMetadata(object):
    def __init__(self, gtf_attrs=None):
        self.chrom = None
        self.start = 0
        self.end = 0
        self.strand = '.'
        self.num_exons = 0
        self.length = 0
        if gtf_attrs is not None:
            for attr in gtf_attrs:
                setattr(self, attr, '')

def get_gtf_metadata(gtf_file, gtf_attrs):
    if gtf_attrs is None:
        gtf_attrs = []
    if 'transcript_id' in gtf_attrs:
        gtf_attrs.remove('transcript_id')
    # read gtf file
    metadata_dict = {}
    for feature in GTFFeature.parse(open(gtf_file)):
        if feature.feature_type != "exon":
            continue
        t_id = feature.attrs["transcript_id"]
        if t_id not in metadata_dict:
            # instantiate new metadata
            m = TranscriptMetadata()
            m.chrom = feature.seqid
            m.strand = feature.strand
            m.start = feature.start
            m.end = feature.end
            for gtf_attr in gtf_attrs:
                setattr(m, gtf_attr, feature.attrs.get(gtf_attr, ''))
            metadata_dict[t_id] = m
        else:
            m = metadata_dict[t_id]
        # update metadata
        m.start = feature.start if feature.start < m.start else m.start
        m.end = feature.end if feature.end > m.end else m.end
        m.length += (feature.end - feature.start)
        m.num_exons += 1
    return metadata_dict

def main():
    # Command line parsing
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    grp = parser.add_mutually_exclusive_group()
    grp.add_argument('--default', action='store_true')
    grp.add_argument('--full', action='store_true')    
    grp.add_argument('-a', '--attr', dest='gtf_attrs', action='append')
    parser.add_argument('gtf_file')
    args = parser.parse_args()
    gtf_file = args.gtf_file
    if not os.path.exists(gtf_file):
        parser.error("gtf file %s not found" % (gtf_file))
    if args.gtf_attrs is not None:
        gtf_attrs = args.gtf_attrs
    elif args.default:
        gtf_attrs = DEFAULT_GTF_ATTRS
    elif args.full:
        gtf_attrs = FULL_GTF_ATTRS
    else:
        gtf_attrs = []
    # collect transcript metadata from GTF file
    logging.info("Reading GTF attributes")
    metadata_dict = get_gtf_metadata(gtf_file, gtf_attrs)
    # output to file
    logging.info("Writing metadata table")
    header_fields = ['transcript_id', 'chrom', 'start', 'end', 'strand', 
                     'num_exons', 'transcript_length'] + gtf_attrs
    print '\t'.join(header_fields)
    for t_id in sorted(metadata_dict):
        m = metadata_dict[t_id]
        fields = [t_id, m.chrom, m.start, m.end, m.strand, m.num_exons, m.length]
        for attr in gtf_attrs:
            fields.append(getattr(m, attr))
        print '\t'.join(map(str,fields))
    return 0

if __name__ == '__main__':
    sys.exit(main())