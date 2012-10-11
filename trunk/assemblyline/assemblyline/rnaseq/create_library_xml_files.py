'''
Created on Aug 7, 2011

@author: mkiyer
'''
import os
import sys
import argparse
import logging
import xml.etree.cElementTree as etree

from libtable import read_library_table_xls
from base import indent_xml

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("library_xls_file")
    parser.add_argument("output_dir")
    args = parser.parse_args()
    libraries = read_library_table_xls(args.library_xls_file)
    for library in libraries.itervalues():
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
    
if __name__ == '__main__': 
    sys.exit(main())
    