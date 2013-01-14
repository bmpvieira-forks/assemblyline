'''
Created on Aug 14, 2011

@author: mkiyer
'''
import argparse
import logging
import subprocess
import sys

# project imports
from assemblyline.rnaseq.lib.inspect import RnaseqLibraryMetrics
from assemblyline.rnaseq.lib.base import detect_read_length

def run_tophat(output_dir, fastq_files, library_metrics_file,
               bowtie_index, num_processors,
               rg_id, rg_sample, rg_library, rg_description, 
               rg_platform_unit, rg_center, rg_platform,
               tophat_args):
    # read library characteristics information
    obj = RnaseqLibraryMetrics.from_file(library_metrics_file)
    library_type = obj.predict_library_type()
    frag_size_mean = int(round(obj.tlen_at_percentile(50.0),0))
    frag_size_stdev = int(round(obj.std(),0))    
    # get fragment size parameters for tophat
    read_length = detect_read_length(fastq_files[0])
    mean_inner_dist = int(frag_size_mean - 2*read_length)
    mate_stdev = frag_size_stdev
    logging.info("Tophat mean_inner_dist=%d and mate_stdev=%d" % 
                 (mean_inner_dist, mate_stdev))
    # setup run
    args = ['tophat',
            "-o", output_dir,
            "--library-type", library_type,
            "-p", num_processors,
            "-r", mean_inner_dist,
            "--mate-std-dev", mate_stdev,
            '--rg-id="%s"' % rg_id,
            '--rg-sample="%s"' % rg_sample,
            '--rg-library="%s"' % rg_library,
            '--rg-description="%s"' % rg_description,
            '--rg-platform-unit="%s"' % rg_platform_unit,
            '--rg-center="%s"' % rg_center,
            '--rg-platform="%s"' % rg_platform]
    for arg in tophat_args:        
        args.extend(arg.split())
    args.append(bowtie_index)
    args.extend(fastq_files)
    logging.debug("\targs: %s" % (map(str, args)))
    retcode = subprocess.call(map(str, args))
    if retcode != 0:
        return retcode

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("--tophat-arg", dest="tophat_args", action="append", default=[])
    parser.add_argument("-p", type=int, dest="num_processors", default=1)
    parser.add_argument("--rg-id", dest="rg_id")
    parser.add_argument("--rg-sample", dest="rg_sample")
    parser.add_argument("--rg-library", dest="rg_library")
    parser.add_argument("--rg-description", dest="rg_description")
    parser.add_argument("--rg-platform-unit", dest="rg_platform_unit")
    parser.add_argument("--rg-center", dest="rg_center")
    parser.add_argument("--rg-platform", dest="rg_platform")
    parser.add_argument("output_dir")
    parser.add_argument("bowtie_index")
    parser.add_argument("library_metrics_file")
    parser.add_argument("fastq_files", nargs="+")
    args = parser.parse_args()
    return run_tophat(output_dir=args.output_dir,
                      fastq_files=args.fastq_files,
                      library_metrics_file=args.library_metrics_file,
                      bowtie_index=args.bowtie_index,
                      num_processors=args.num_processors, 
                      rg_id=args.rg_id,
                      rg_sample=args.rg_sample, 
                      rg_library=args.rg_library, 
                      rg_description=args.rg_description, 
                      rg_platform_unit=args.rg_platform_unit, 
                      rg_center=args.rg_center, 
                      rg_platform=args.rg_platform,
                      tophat_args=args.tophat_args)
    
if __name__ == '__main__':
    sys.exit(main())