'''
Created on Dec 2, 2012

@author: mkiyer
'''
import logging
import argparse
import sys
import os

import assemblyline.rnaseq.lib.config as config
from assemblyline.rnaseq.lib.libtable import Library

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("output_dir")
    args = parser.parse_args()
    library_xml_file = os.path.join(args.output_dir, config.LIBRARY_XML_FILE)
    if not os.path.exists(library_xml_file):
        logging.error("Library xml file not found")
        return config.JOB_ERROR
    library = Library.from_xml_file(library_xml_file)
    if not library.is_valid():
        logging.error("Library not valid")
        return config.JOB_ERROR
    results = config.RnaseqResults(library, args.output_dir)
    valid, missing_files = results.validate()
    if not valid:
        valid = (len(missing_files) == 1 and 
                 missing_files[0] == results.job_done_file)
    if valid:
        open(results.job_done_file, "w").close()
    return config.JOB_SUCCESS
    
if __name__ == '__main__':
    sys.exit(main())