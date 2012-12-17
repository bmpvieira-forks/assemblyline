'''
Created on Oct 15, 2012

@author: mkiyer
'''
import logging
import argparse
import sys
import os

import assemblyline.rnaseq.lib.config as config
from assemblyline.rnaseq.lib.libtable import Library, read_library_table_xls

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("library_xls_file")
    parser.add_argument("root_dir")
    args = parser.parse_args()
    # read library file
    logging.info("Reading library table '%s'" % (args.library_xls_file))
    libraries = read_library_table_xls(args.library_xls_file)
    # build list of all parameters
    params = set()
    for library in libraries.itervalues():
        params.update(library.params.keys())
    sorted_params = sorted(params)
    # output table
    logging.info("Generating library table")
    header_fields = []
    header_fields.extend(Library.fields)
    header_fields.extend(sorted_params)
    header_fields.extend(["bam_file", "gtf_file"])
    print '\t'.join(header_fields)
    for library in libraries.itervalues():
        output_dir = os.path.join(args.root_dir, library.library_id)
        results = config.RnaseqResults(library, output_dir)
        if not results.validate():
            setattr(library, 'description', 'invalid')
        else:
            setattr(library, 'description', 'valid')
        library.read1_files = ','.join(library.read1_files)
        library.read2_files = ','.join(library.read2_files)
        fields = []
        for field_name in Library.fields:
            fields.append(getattr(library, field_name))
        for param in sorted_params:
            fields.append(library.params.get(param, "na"))
        fields.append(results.tophat_bam_file)
        fields.append(results.cufflinks_gtf_file)
        print '\t'.join(fields)
    return config.JOB_SUCCESS
    
if __name__ == '__main__':
    sys.exit(main())