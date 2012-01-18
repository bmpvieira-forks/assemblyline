'''
Created on Nov 28, 2011

@author: mkiyer
'''
import logging
import argparse
import os
import sys

from assemblyline.lib.sampletable import LibInfo
from assemblyline.lib import gtf

def main():
    #
    # Command line parsing
    # 
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', dest="output_file", default=None)
    parser.add_argument('lib_table')
    args = parser.parse_args()
    # open output file
    if args.output_file is not None:
        outfh = open(args.output_file, "w")
    else:
        outfh = sys.stdout
    # parse library table
    logging.info("Parsing library table")
    for libinfo in LibInfo.from_file(args.lib_table):
        # exclude invalid libraries
        if (not libinfo.is_valid()) or (not libinfo.for_assembly):
            logging.warning("\tskipping cohort=%s patient=%s sample=%s lib=%s lanes=%s" % 
                            (libinfo.cohort, libinfo.patient, libinfo.sample, 
                             libinfo.library, libinfo.lanes))
            continue
        logging.info("\tadding cohort=%s patient=%s sample=%s lib=%s lanes=%s" % 
                     (libinfo.cohort, libinfo.patient, libinfo.sample, 
                      libinfo.library, libinfo.lanes))        
        # parse and annotate GTF files
        for feature in gtf.GTFFeature.parse(open(libinfo.cufflinks_gtf_file)):
            feature.attrs['cohort'] = libinfo.cohort
            feature.attrs['sample'] = libinfo.sample
            feature.attrs['library'] = libinfo.library
            print >>outfh, str(feature)
    # cleanup
    if args.output_file is not None:
        outfh.close()

if __name__ == '__main__':
    sys.exit(main())
