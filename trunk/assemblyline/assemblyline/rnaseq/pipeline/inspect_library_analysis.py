'''
Created on Dec 5, 2012

@author: mkiyer
'''
import argparse
import logging
import sys

# local import
import pysam

# project imports
import assemblyline.rnaseq.lib.config as config
from assemblyline.rnaseq.lib.libtable import FRAGMENT_LAYOUT_PAIRED
from assemblyline.rnaseq.lib.inspect import inspect_sr_sam, inspect_pe_sam, RnaseqLibraryMetrics

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("-n", dest="n", type=int, default=1000000)
    parser.add_argument("--strand-spec-frac", dest="strand_spec_frac",
                        type=float, default=0.90)
    parser.add_argument("--fragment-layout", dest="fragment_layout", 
                        default=FRAGMENT_LAYOUT_PAIRED)
    parser.add_argument("--min-frag-size", dest="min_frag_size", 
                        type=int, default=0)
    parser.add_argument("--max-frag-size", dest="max_frag_size", 
                        type=int, default=1000)
    parser.add_argument("--frag-size-mean", dest="frag_size_mean_default", 
                        type=float, default=200)
    parser.add_argument("--frag-size-stdev", dest="frag_size_stdev_default", 
                        type=float, default=20)
    parser.add_argument("sam_file")
    parser.add_argument("output_file")
    args = parser.parse_args()
    samfh = pysam.Samfile(args.sam_file, "r")
    obj = RnaseqLibraryMetrics(args.min_frag_size, 
                               args.max_frag_size,
                               args.strand_spec_frac)
    if args.fragment_layout == FRAGMENT_LAYOUT_PAIRED:
        obj.from_stream(inspect_pe_sam(samfh, args.n))
    else:
        obj.from_stream(inspect_sr_sam(samfh, args.n))
    samfh.close()
    # if the number of data points is too few, use a default distribution
    num_frag_size_samples = obj.num_frag_size_samples()
    num_strand_samples = obj.read1_count
    if (num_frag_size_samples < config.MIN_INSPECT_SAMPLES):
        logging.warning("Found %d (<%d) samples found during profiling, "
                        "so cannot reliably estimate the fragment "
                        "size distribution.  Using defaults instead." % 
                        (num_frag_size_samples, config.MIN_INSPECT_SAMPLES))
        obj.clear_frag_size()
        obj.frag_size_from_random(args.frag_size_mean_default,
                                  args.frag_size_stdev_default,
                                  args.n)
    if (num_strand_samples < config.MIN_INSPECT_SAMPLES):
        logging.warning("Found %d (<%d) samples found during profiling, "
                        "so cannot reliably estimate strandedness. "
                        "Using defaults instead." % 
                        (num_strand_samples, config.MIN_INSPECT_SAMPLES))
        obj.clear_strandedness()
        obj.strandedness_from_random(args.n)
    # write the distribution
    obj.to_file(open(args.output_file, "w"))
    return config.JOB_SUCCESS

if __name__ == '__main__': 
    sys.exit(main())