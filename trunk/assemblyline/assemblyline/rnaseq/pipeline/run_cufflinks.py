'''
Created on Aug 14, 2011

@author: mkiyer
'''
import argparse
import logging
import subprocess
import sys
import os

from assemblyline.rnaseq.lib.inspect import RnaseqLibraryMetrics
import assemblyline.rnaseq.lib.config as config

def run_cufflinks(bam_file,
                  output_dir,
                  output_label,
                  library_type,
                  cufflinks_args,
                  num_processors=1,
                  learn_frag_size_dist=True,
                  library_metrics_file=None):
    # create output dir if necessary
    if not os.path.exists(output_dir):
        logging.info("Creating directory: %s" % (output_dir))
        os.makedirs(output_dir)
    # read rnaseq characteristics info
    obj = RnaseqLibraryMetrics.from_file(open(library_metrics_file))
    # predict library type
    predicted_library_type = obj.predict_library_type(config.STRAND_SPECIFIC_CUTOFF_FRAC)
    if library_type != predicted_library_type:
        logging.warning("Library type %s does not agree with prediction %s" % 
                        (library_type, predicted_library_type))
    args = ['cufflinks',
            "-o", output_dir,
            "-p", num_processors,
            "--library-type", library_type,
            "-L", output_label]
    for arg in cufflinks_args:        
        args.extend(arg.split())
    if not learn_frag_size_dist:
        # get fragment size parameters
        args.extend(["-m", int(obj.tlen_at_percentile(50.0)),
                     "-s", int(round(obj.std(),0))])
    args.append(bam_file)
    logging.debug("Cufflinks command line: %s" % (" ".join(map(str, args))))
    retcode = subprocess.call(map(str, args))
    return retcode

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", type=int, dest="num_processors", default=1)
    parser.add_argument("-L", dest="output_label", default="CUFF")
    parser.add_argument("--library-type", dest="library_type", default="fr-unstranded")
    parser.add_argument("--learn-frag-size", dest="learn_frag_size", action="store_true", default=False)
    parser.add_argument("--cufflinks-arg", dest="cufflinks_args", action="append", default=[])
    parser.add_argument("bam_file")
    parser.add_argument("output_dir")
    parser.add_argument("library_metrics_file")
    args = parser.parse_args()
    return run_cufflinks(bam_file=args.bam_file,
                         output_dir=args.output_dir,
                         output_label=args.output_label,
                         library_type=args.library_type,
                         cufflinks_args=args.cufflinks_args,
                         num_processors=args.num_processors,
                         learn_frag_size_dist=args.learn_frag_size,
                         library_metrics_file=args.library_metrics_file)
    
if __name__ == '__main__':
    sys.exit(main())
