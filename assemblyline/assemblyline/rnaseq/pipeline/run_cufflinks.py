'''
Created on Aug 14, 2011

@author: mkiyer
'''
import argparse
import logging
import subprocess
import sys

from assemblyline.rnaseq.lib.fragment_size_distribution import FragmentSizeDistribution

def run_cufflinks(bam_file,
                  output_dir,
                  output_label,
                  library_type,
                  cufflinks_args,
                  num_processors=1,
                  learn_frag_size_dist=True,
                  frag_size_dist_file=None,
                  cufflinks_bin="cufflinks"):
    args = [cufflinks_bin,
            "-o", output_dir,
            "-p", num_processors,
            "--library-type", library_type,
            "-L", output_label]
    for arg in cufflinks_args:        
        args.extend(arg.split())
    if not learn_frag_size_dist:
        # get fragment size parameters
        frag_size_dist = FragmentSizeDistribution.from_file(open(frag_size_dist_file))
        args.extend(["-m", int(frag_size_dist.isize_at_percentile(50.0)),
                     "-s", int(round(frag_size_dist.std(),0))])
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
    parser.add_argument("--cufflinks-bin", dest="cufflinks_bin", default="cufflinks")
    parser.add_argument("--cufflinks-arg", dest="cufflinks_args", action="append", default=[])
    parser.add_argument("bam_file")
    parser.add_argument("output_dir")
    parser.add_argument("frag_size_dist_file")
    args = parser.parse_args()
    return run_cufflinks(bam_file=args.bam_file,
                         output_dir=args.output_dir,
                         output_label=args.output_label,
                         library_type=args.library_type,
                         cufflinks_args=args.cufflinks_args,
                         num_processors=args.num_processors,
                         learn_frag_size_dist=args.learn_frag_size,
                         frag_size_dist_file=args.frag_size_dist_file,
                         cufflinks_bin=args.cufflinks_bin)

    
if __name__ == '__main__':
    sys.exit(main())
