'''
Created on Dec 17, 2013

@author: mkiyer
'''
import sys
import argparse
import logging

from assemblyline.utils.gtf_merge import Feature

def genepred_to_gtf(bed_file, source=None):
    # read all features
    attr_fields=['gene_name']    
    for f in Feature.parse_genepred(open(bed_file), attr_fields):
        for gtf_feature in f.to_gtf_features(source=source):
            yield str(gtf_feature)

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument('--source', dest='source', default=None)
    parser.add_argument('genepred_file')
    args = parser.parse_args()
    for line in genepred_to_gtf(args.genepred_file, args.source):
        print line

if __name__ == '__main__':
    sys.exit(main())