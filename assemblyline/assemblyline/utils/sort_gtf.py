'''
Created on Nov 4, 2010

@author: mkiyer
'''
import os
import sys
import logging
import argparse

import assemblyline
from assemblyline.lib.gtf import sort_gtf

def main():
    # setup logging
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    logging.info("AssemblyLine %s" % (assemblyline.__version__))
    logging.info("----------------------------------")
    # parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("input_gtf_file")
    parser.add_argument("output_gtf_file")
    args = parser.parse_args()
    # check command line
    if not os.path.exists(args.input_gtf_file):
        parser.error("input gtf file %s not found" % (args.input_gtf_file))
    logging.info("Parameters:")
    logging.info("input gtf file:  %s" % (args.input_gtf_file))
    logging.info("output gtf file: %s" % (args.output_gtf_file))
    logging.info("Sorting")
    sort_gtf(args.input_gtf_file, args.output_gtf_file)    
    logging.info("Done")
    return 0

if __name__ == '__main__':
    sys.exit(main())
