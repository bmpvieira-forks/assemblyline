'''
Created on Dec 10, 2012

@author: mkiyer
'''
import os
import sys
import argparse
import logging
import subprocess

from assemblyline.rnaseq.lib.inspect import RnaseqLibraryMetrics

def get_picard_strand_specificity(library_type):
    if library_type == "fr-firststrand":
        return "SECOND_READ_TRANSCRIPTION_STRAND"
    elif library_type == "fr-secondstrand":        
        return "FIRST_READ_TRANSCRIPTION_STRAND"
    else:
        return "NONE"

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("--tmp-dir", dest="tmp_dir", default="/tmp")
    parser.add_argument("--picard-dir", dest="picard_dir", default="")
    parser.add_argument("tophat_bam_file")
    parser.add_argument("library_metrics_file")
    parser.add_argument("refflat_file")
    parser.add_argument("ribosomal_intervals_file")
    parser.add_argument("reference_fasta")
    parser.add_argument("output_file")
    parser.add_argument("chart_file")
    args = parser.parse_args()    
    logging.debug("Predicting library type from metrics file")
    obj = RnaseqLibraryMetrics.from_file(args.library_metrics_file)
    predicted_library_type = obj.predict_library_type()
    logging.debug("Predicted library type: %s" % (predicted_library_type))
    picard_args = ["java", "-Xmx4g", "-jar", 
                   os.path.join(args.picard_dir, "CollectRnaSeqMetrics.jar"),
                   "INPUT=%s" % (args.tophat_bam_file),
                   "REF_FLAT=%s" % (args.refflat_file),
                   "RIBOSOMAL_INTERVALS=%s" % (args.ribosomal_intervals_file),
                   "STRAND_SPECIFICITY=%s" % get_picard_strand_specificity(predicted_library_type),
                   "REFERENCE_SEQUENCE=%s" % (args.reference_fasta),
                   "OUTPUT=%s" % (args.output_file),
                   "CHART_OUTPUT=%s" % (args.chart_file),
                   "TMP_DIR=%s" % args.tmp_dir]
    logging.debug("\targs: %s" % (' '.join(map(str, picard_args))))
    subprocess.call(picard_args)

if __name__ == '__main__':
    sys.exit(main())