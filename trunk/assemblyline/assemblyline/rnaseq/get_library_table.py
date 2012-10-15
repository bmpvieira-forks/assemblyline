'''
Created on Oct 15, 2012

@author: mkiyer
'''
import logging
import argparse
import sys

import assemblyline.rnaseq.lib.config as config
from assemblyline.rnaseq.lib.libtable import Library, read_library_table_xls

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("config_xml_file")
    parser.add_argument("server_name")
    parser.add_argument("library_xls_file")
    args = parser.parse_args()
    # read configuration file
    logging.info("Reading pipeline configuration file '%s'" % (args.config_xml_file))
    pipeline = config.PipelineConfig.from_xml(args.config_xml_file)
    # get server information
    if args.server_name not in pipeline.servers:
        logging.error("Server %s not found" % (args.server_name))
    server = pipeline.servers[args.server_name]
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
        results = config.RnaseqResults(library, server.output_dir)
        if not results.validate():
            continue
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