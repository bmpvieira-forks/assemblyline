'''
Created on Dec 5, 2012

@author: mkiyer
'''
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
from assemblyline.rnaseq.lib.inspect import RnaseqLibraryMetrics
from assemblyline.rnaseq.lib.libtable import FRAGMENT_LAYOUT_PAIRED, FRAGMENT_LAYOUT_SINGLE
import assemblyline.rnaseq.lib.config as config

import assemblyline.rnaseq.pipeline
_pipeline_dir = assemblyline.rnaseq.pipeline.__path__[0]

def inspect_rnaseq_library(bowtie_index,
                           result_file,
                           frag_size_dist_plot_file,
                           fastq_files,
                           strand_spec_frac,
                           frag_size_mean_default,
                           frag_size_stdev_default,
                           max_num_samples,
                           min_frag_size,
                           max_frag_size,
                           num_processors):
    bowtie_threads = max(1, num_processors-1)
    args = ["bowtie", "-q", "-S",
            "-p", str(bowtie_threads),
            "-v", 2, "-k", 1, "-m", 1]
    # different approaches for single-read or paired-end
    if len(fastq_files) == 1:
        # for single read we cannot readily estimate the fragment size 
        # distribution, but can get strandedness
        fragment_layout = FRAGMENT_LAYOUT_SINGLE
        args.extend([bowtie_index, fastq_files[0]])
    else:
        # setup bowtie process for paired-end mapping
        fragment_layout = FRAGMENT_LAYOUT_PAIRED
        args.extend(["--minins", min_frag_size,
                     "--maxins", max_frag_size,
                     bowtie_index,
                     "-1", fastq_files[0],
                     "-2", fastq_files[1]])
    logging.debug("bowtie args: %s" % (' '.join(map(str, args))))
    aln_p = subprocess.Popen(map(str,args), stdout=subprocess.PIPE)
    # setup inspection process
    args = ["python", 
            os.path.join(_pipeline_dir, "inspect_library_analysis.py"),
            "--strand-spec-frac", strand_spec_frac,
            "--fragment-layout", fragment_layout,
            "-n", max_num_samples,
            "--min-frag-size", min_frag_size,
            "--max-frag-size", max_frag_size,
            "--frag-size-mean", frag_size_mean_default,
            "--frag-size-stdev", frag_size_stdev_default,
            "-",
            result_file]
    logging.debug("inspect analysis args: %s" % (' '.join(map(str, args))))
    retcode = subprocess.call(map(str,args), stdin=aln_p.stdout)    
    aln_p.terminate()
    if retcode != 0:
        logging.error("inspect analysis failed")
        # cleanup output file
        if os.path.exists(result_file):
            os.remove(result_file)
        return config.JOB_ERROR
    # plot fragment size distribution
    obj = RnaseqLibraryMetrics.from_file(result_file)
    logging.info("Fragment size samples=%d mean=%f std=%f median=%d mode=%d" % 
                 (obj.num_frag_size_samples(), obj.mean(), obj.std(), 
                  obj.tlen_at_percentile(50.0), obj.mode()))    
    logging.info("Strandedness frac=%f samples=%d rev=%d predicted library type=%s" %
                 (obj.read1_rev_count / float(obj.read1_count), 
                  obj.read1_count, obj.read1_rev_count, obj.predict_library_type()))
    logging.info("Writing metrics file")
    obj.to_file(open(result_file, "w"))
    logging.info("Plotting fragment size distribution")
    obj.plot_frag_size_dist(frag_size_dist_plot_file, max_percentile=99.9)
    logging.info("Done")
    return config.JOB_SUCCESS

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("--strand-spec-frac", dest="strand_spec_frac",
                        type=float, default=0.90)
    parser.add_argument("--min-frag-size", dest="min_frag_size", 
                        type=int, default=0)
    parser.add_argument("--max-frag-size", dest="max_frag_size", 
                        type=int, default=1000)
    parser.add_argument("--frag-size-mean", dest="frag_size_mean_default", 
                        type=float, default=200)
    parser.add_argument("--frag-size-stdev", dest="frag_size_stdev_default", 
                        type=float, default=20)
    parser.add_argument("-n", dest="max_num_samples", type=int, default=1000000)
    parser.add_argument("-p", dest="num_processors", type=int, default=1)
    parser.add_argument("bowtie_index")
    parser.add_argument("result_file")
    parser.add_argument("frag_size_dist_plot_file")
    parser.add_argument("fastq_files", nargs="+")
    args = parser.parse_args()
    return inspect_rnaseq_library(args.bowtie_index,
                                  args.result_file,
                                  args.frag_size_dist_plot_file,
                                  args.fastq_files,
                                  strand_spec_frac=args.strand_spec_frac,
                                  frag_size_mean_default=args.frag_size_mean_default,
                                  frag_size_stdev_default=args.frag_size_stdev_default,
                                  max_num_samples=args.max_num_samples,
                                  min_frag_size=args.min_frag_size,
                                  max_frag_size=args.max_frag_size,
                                  num_processors=args.num_processors)

if __name__ == '__main__': 
    sys.exit(main())