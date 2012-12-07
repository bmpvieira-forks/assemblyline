'''
Created on Aug 10, 2011

@author: mkiyer
'''
import argparse
import logging
import sys
import os
import subprocess

# project imports
from assemblyline.rnaseq.lib.fragment_size_distribution import FragmentSizeDistribution
import assemblyline.rnaseq.lib.config as config
import assemblyline.rnaseq.pipeline
_pipeline_dir = assemblyline.rnaseq.pipeline.__path__[0]

def estimate_fragment_size_distribution(bowtie_index, 
                                        output_file,
                                        plot_file,
                                        fastq_files,
                                        frag_size_mean_default=200,
                                        frag_size_stdev_default=20,
                                        num_samples=1000000,
                                        min_fragment_size=0,
                                        max_fragment_size=1000,
                                        trim5=0,
                                        num_processors=1,
                                        bowtie_bin="bowtie"):
    # for single read we cannot easily estimate the fragment size 
    # distribution so use defaults
    if len(fastq_files) == 1:
        frag_size_dist = FragmentSizeDistribution.from_random(frag_size_mean_default, 
                                                              frag_size_stdev_default,
                                                              min_fragment_size, 
                                                              max_fragment_size, 
                                                              1000000)
        frag_size_dist.to_file(open(output_file, "w"))
    else:
        # setup bowtie process
        bowtie_threads = max(1, num_processors-1)
        args = [bowtie_bin, "-q", "-S",
                "-p", str(bowtie_threads),
                "-v", 2, "-k", 1, "-m", 1, 
                "--minins", min_fragment_size,
                "--maxins", max_fragment_size,
                bowtie_index,
                "-1", fastq_files[0],
                "-2", fastq_files[1]]
        logging.debug("Fragment size bowtie args: %s" % (' '.join(map(str, args))))
        aln_p = subprocess.Popen(map(str,args), stdout=subprocess.PIPE)
        # setup fragment size estimation process
        args = [sys.executable, os.path.join(_pipeline_dir, "process_fragment_size_sam.py"),
                "--min-fragment-size", min_fragment_size,
                "--max-fragment-size", max_fragment_size,
                "--trim5", trim5,
                "-n", num_samples,
                "-", fastq_files[0], output_file]
        logging.debug("Fragment size processing args: %s" % (' '.join(map(str, args))))
        retcode = subprocess.call(map(str,args), stdin=aln_p.stdout)    
        aln_p.terminate()
        if retcode != 0:
            logging.error("Fragment size processing script failed")
            # cleanup output file
            if os.path.exists(output_file):
                os.remove(output_file)
            return config.JOB_ERROR
        # read the distribution
        frag_size_dist = FragmentSizeDistribution.from_file(open(output_file, "r"))
        # if the number of data points is too few, use a default distribution
        if frag_size_dist.n < config.MIN_FRAG_SIZE_SAMPLES:
            logging.warning("Found %d (<%d) samples found during profiling, "
                            "so cannot reliably estimate the fragment "
                            "size distribution.  Using defaults instead." % 
                            (frag_size_dist.n, config.MIN_FRAG_SIZE_SAMPLES))
            frag_size_dist = FragmentSizeDistribution.from_random(frag_size_mean_default, 
                                                                  frag_size_stdev_default,
                                                                  min_fragment_size, 
                                                                  max_fragment_size, 
                                                                  1000000)
            frag_size_dist.to_file(open(output_file, "w"))
    # plot the distribution
    frag_size_dist.plot(plot_file, max_percentile=99.9)
    return config.JOB_SUCCESS    

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("--bowtie-bin", dest="bowtie_bin",
                        default="bowtie")
    parser.add_argument("--min-fragment-size", dest="min_fragment_size", 
                        type=int, default=0)
    parser.add_argument("--max-fragment-size", dest="max_fragment_size", 
                        type=int, default=1000)
    parser.add_argument("--trim5", dest="trim5", 
                        type=int, default=0)
    parser.add_argument("--default-mean", dest="frag_size_mean_default", 
                        type=float, default=200)
    parser.add_argument("--default-stdev", dest="frag_size_stdev_default", 
                        type=float, default=20)
    parser.add_argument("-n", dest="num_samples", type=int, default=1000000)
    parser.add_argument("-p", dest="num_processors", type=int, default=1)
    parser.add_argument("bowtie_index")
    parser.add_argument("output_file")
    parser.add_argument("plot_file")
    parser.add_argument("fastq_files", nargs="+")
    args = parser.parse_args()
    return estimate_fragment_size_distribution(args.bowtie_index,
                                               args.output_file,
                                               args.plot_file,
                                               args.fastq_files,
                                               frag_size_mean_default=args.frag_size_mean_default,
                                               frag_size_stdev_default=args.frag_size_stdev_default,
                                               num_samples=args.num_samples,
                                               min_fragment_size=args.min_fragment_size,
                                               max_fragment_size=args.max_fragment_size,
                                               trim5=args.trim5,
                                               num_processors=args.num_processors,
                                               bowtie_bin=args.bowtie_bin)

if __name__ == '__main__': 
    sys.exit(main())