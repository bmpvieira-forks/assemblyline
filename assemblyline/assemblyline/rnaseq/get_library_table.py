'''
Created on Oct 15, 2012

@author: mkiyer
'''
import logging
import argparse
import sys
import os

import assemblyline.rnaseq.lib.config as config
import assemblyline.rnaseq.lib.picard as picard
from assemblyline.rnaseq.lib.libtable import Library, read_library_table_xls
from assemblyline.rnaseq.lib.base import get_fastqc_total_sequences

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
    header_fields.extend(["total_fragments",
                          "total_aligned_reads", 
                          "mean_read_length", 
                          'aligned_bases',
                          'pct_mrna_bases',
                          'pct_intronic_bases',
                          'pct_intergenic_bases',
                          'median_5prime_to_3prime_bias',
                          'median_cv_coverage',
                          'bam_file', 
                          'gtf_file'])
    print '\t'.join(header_fields)
    for library_id in sorted(libraries):
        library = libraries[library_id]
        output_dir = os.path.join(args.root_dir, library.library_id)
        config_xml_file = os.path.join(output_dir, config.CONFIG_XML_FILE)
        if not os.path.exists(config_xml_file):
            logging.error("Configuration xml file not found")
            continue
        pipeline = config.PipelineConfig.from_xml(args.config_xml_file)
        results = config.RnaseqResults(library, pipeline, output_dir)
        is_valid, missing_files = results.validate()
        if not is_valid:
            setattr(library, 'description', 'invalid')
        else:
            setattr(library, 'description', 'valid')
        library.read1_files = ','.join(library.read1_files)
        library.read2_files = ','.join(library.read2_files)
        library.bam_files = ','.join(library.bam_files)
        fields = []
        for field_name in Library.fields:
            fields.append(getattr(library, field_name))
        for param in sorted_params:
            fields.append(library.params.get(param, "na"))
        # results
        if os.path.exists(results.fastqc_data_files[0]):
            total_frags = get_fastqc_total_sequences(results.fastqc_data_files[0])
            fields.append(total_frags)
        else:
            fields.append('na')
        if os.path.exists(results.alignment_summary_metrics):
            obj = picard.AlignmentSummaryMetrics(results.alignment_summary_metrics)
            fields.append(obj.get_total_reads())
            fields.append(obj.get_mean_read_length()) 
        else:
            fields.extend(['na', 'na'])
        if os.path.exists(results.rnaseq_metrics):
            metrics_dict = picard.get_rnaseq_metrics(results.rnaseq_metrics)
            fields.extend([metrics_dict['PF_ALIGNED_BASES'],
                           metrics_dict['PCT_MRNA_BASES'],
                           metrics_dict['PCT_INTRONIC_BASES'],
                           metrics_dict['PCT_INTERGENIC_BASES'],
                           metrics_dict['MEDIAN_5PRIME_TO_3PRIME_BIAS'],
                           metrics_dict['MEDIAN_CV_COVERAGE']])
        if os.path.exists(results.tophat_bam_file):
            fields.append(results.tophat_bam_file)
        else:
            fields.append('na')
        if os.path.exists(results.cufflinks_gtf_file):
            fields.append(results.cufflinks_gtf_file)
        else:
            fields.append('na')
        print '\t'.join(map(str, fields))
    return config.JOB_SUCCESS
    
if __name__ == '__main__':
    sys.exit(main())