'''
Created on Aug 10, 2011

@author: mkiyer
'''
import argparse
import logging
import sys

# local import
import pysam

# project imports
import assemblyline.rnaseq.config as config
from assemblyline.rnaseq.fragment_size_distribution import FragmentSizeDistribution

def process_fragment_size_sam(sam_file,
                              fastq_file,
                              output_file,
                              num_samples=1000000,
                              min_fragment_size=0,
                              max_fragment_size=1000,
                              trim5=0):
    # TODO: incorporate trim5
    samfh = pysam.Samfile(sam_file, "r")
    d = FragmentSizeDistribution.from_sam(samfh, 
                                          min_isize=min_fragment_size,
                                          max_isize=max_fragment_size,
                                          max_samples=num_samples)
    logging.info("Fragment size samples=%d mean=%f std=%f median=%d mode=%d" % 
                 (d.n, d.mean(), d.std(), d.isize_at_percentile(50.0), d.mode()))    
    d.to_file(open(output_file, "w"))
    return config.JOB_SUCCESS

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("--min-fragment-size", dest="min_fragment_size", 
                        type=int, default=0)
    parser.add_argument("--max-fragment-size", dest="max_fragment_size", 
                        type=int, default=1000)
    parser.add_argument("--trim5", dest="trim5", 
                        type=int, default=0)
    parser.add_argument("-n", dest="num_samples", type=int, default=1000000)
    parser.add_argument("sam_file")
    parser.add_argument("fastq_file")
    parser.add_argument("output_file")    
    args = parser.parse_args()
    return process_fragment_size_sam(args.sam_file,
                                     args.fastq_file,
                                     args.output_file,
                                     num_samples=args.num_samples,
                                     min_fragment_size=args.min_fragment_size,
                                     max_fragment_size=args.max_fragment_size,
                                     trim5=args.trim5)

if __name__ == '__main__': 
    sys.exit(main())