'''
Created on Dec 19, 2012

@author: mkiyer
'''
import sys
import os
import logging
import argparse

# projects
import assemblyline.rnaseq.lib.config as config
from assemblyline.rnaseq.lib.libtable import Library
from assemblyline.rnaseq.lib.inspect import RnaseqLibraryMetrics
from assemblyline.rnaseq.lib.libtable import FR_UNSTRANDED

import assemblyline.rnaseq.pipeline
_pipeline_dir = assemblyline.rnaseq.pipeline.__path__[0]

STRAND_COLOR_MAP = {".": "0,128,0",
                    "+": "255,0,0",
                    "-": "0,0,255"}

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("--baseurl", dest="baseurl")
    parser.add_argument("library_dir")
    args = parser.parse_args()
    if not os.path.exists(args.library_dir):
        parser.error("library_dir '%s' not found" % (args.library_dir))
    if not os.path.isdir(args.library_dir):
        parser.error("library_dir '%s' not a directory" % (args.library_dir))        
    library_dir = os.path.abspath(args.library_dir)
    # read library file
    library_xml_file = os.path.join(library_dir, config.LIBRARY_XML_FILE)
    library = Library.from_xml_file(library_xml_file)
    # get results
    results = config.RnaseqResults(library, library_dir)
    # predict library type
    obj = RnaseqLibraryMetrics.from_file(open(results.library_metrics_file))
    library_type = obj.predict_library_type(config.STRAND_SPECIFIC_CUTOFF_FRAC)
    #
    # coverage files (bigwig)
    # 
    if library_type == FR_UNSTRANDED:
        strand_files = [(".", results.coverage_bigwig_prefix + ".bw")]
    else:
        strand_files = [("+", results.coverage_bigwig_prefix + "_pos.bw"),
                        ("-", results.coverage_bigwig_prefix + "_neg.bw")] 
    for strand,bigwig_file in strand_files:
        ucsc_url = "%s%s" % (args.baseurl, bigwig_file)
        track_name = "cov_%s_%s" % (results.library_id, strand)
        track_desc = "Coverage for %s strand=[%s] (RPM)" % (results.library_id, strand)
        track_line = ' '.join(['track type=bigWig',
                               'name="%s"' % (track_name),
                               'description="%s"' % (track_desc),
                               'visibility=full',
                               'color=%s' % (STRAND_COLOR_MAP[strand]),
                               'autoScale=on',
                               'maxHeightPixels=64:64:11',
                               'bigDataUrl="%s"' % (ucsc_url)])
        print track_line
    #
    # junction bigbed file
    #
    ucsc_url = "%s%s" % (args.baseurl, results.junctions_bigbed_file)
    track_name = "junc_%s" % (results.library_id)
    track_desc = "Splice junctions for %s" % (results.library_id)
    track_line = ' '.join(['track type=bigBed',
                           'name="%s"' % (track_name),
                           'description="%s"' % (track_desc),
                           'visibility=pack',
                           'bigDataUrl="%s"' % (ucsc_url)])
    print track_line
    #
    # variant calls (VCF)
    #
    ucsc_url = "%s%s" % (args.baseurl, results.varscan_snv_bgzip_file)
    track_name = "snv_%s" % (results.library_id)
    track_desc = "SNV for %s" % (results.library_id)
    track_line = ' '.join(['track type=vcfTabix',
                           'name="%s"' % (track_name),
                           'description="%s"' % (track_desc),
                           'visibility=pack',
                           'bigDataUrl="%s"' % (ucsc_url)])
    print track_line

if __name__ == '__main__':
    sys.exit(main())