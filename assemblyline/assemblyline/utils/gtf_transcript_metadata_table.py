'''
Created on Feb 24, 2013

@author: mkiyer
'''
import logging
import argparse
import os
import sys

from make_expression_matrix import get_gtf_metadata

def main():
    # Command line parsing
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument('--group-by', dest="group_by",
                        choices=["gene_id", "transcript_id"],
                        default="gene_id")
    parser.add_argument('--omit', dest="omit_list",
                        default="exon_number",
                        help="default=%(default)s")
    parser.add_argument('gtf_file')
    args = parser.parse_args()
    if not os.path.exists(args.gtf_file):
        parser.error("gtf file %s not found" % (args.gtf_file))
    omit_list = args.omit_list.split(",")
    # get gene metadata
    logging.info("Reading metadata")
    metadata_fields, metadata_dict = \
        get_gtf_metadata(args.gtf_file, omit_attrs=omit_list,
                         group_by=args.group_by)
    # write matrix file
    logging.info("Writing table")
    header_fields = list(metadata_fields)
    print '\t'.join(header_fields)
    for k in sorted(metadata_dict):
        fields = list(map(str, metadata_dict[k]))
        print '\t'.join(fields)
    return 0

if __name__ == '__main__':
    sys.exit(main())