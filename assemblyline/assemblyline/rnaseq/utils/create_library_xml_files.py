'''
Created on Aug 7, 2011

@author: mkiyer
'''
import os
import sys
import argparse
import logging
import xml.etree.cElementTree as etree

import assemblyline.rnaseq.lib.config as config 
from assemblyline.rnaseq.lib.base import indent_xml
from assemblyline.rnaseq.lib.libtable import read_library_table_xls

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("config_xml_file")
    parser.add_argument("server_name")    
    parser.add_argument("library_xls_file")
    parser.add_argument("output_dir")
    args = parser.parse_args()
    # read configuration file
    logging.info("Reading pipeline configuration file '%s'" % (args.config_xml_file))
    pipeline = config.PipelineConfig.from_xml(args.config_xml_file)
    if not pipeline.is_valid(args.server_name):
        logging.error("Pipeline config not valid")
        return config.JOB_ERROR
    server = pipeline.servers[args.server_name]
    # read library file
    logging.info("Reading library table '%s'" % (args.library_xls_file))
    libraries = read_library_table_xls(args.library_xls_file)
    logging.info("Generating XML files at output directory '%s'" % (args.output_dir))
    for library in libraries.itervalues():
        # search for library sequence files
        success = config.resolve_library_sequence_files(server, library)
        if not success:
            logging.error("Library %s sequence files not found" % (library.library_id))
            continue
        # check for valid data structure
        if not library.is_valid():
            logging.error("Library %s not valid" % (library.library_id))
            continue
        # add patient        
        root = etree.Element("library")
        library.to_xml(root)
        # indent for pretty printing
        indent_xml(root)
        # write to file
        output_file = os.path.join(args.output_dir, "%s.xml" % (library.library_id))        
        f = open(output_file, "w")
        print >>f, etree.tostring(root)
        f.close()  
    logging.info("Done")
    
if __name__ == '__main__': 
    sys.exit(main())
    