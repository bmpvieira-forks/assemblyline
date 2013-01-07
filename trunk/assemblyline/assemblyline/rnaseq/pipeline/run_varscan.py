'''
Created on Feb 1, 2012

@author: mkiyer
'''
import logging
import argparse
import subprocess
import sys
import os

def run_varscan(ref_fa, 
                bam_file, 
                output_file,
                varscan_jar,
                varscan_func,
                varscan_args):
    # run varscan snv
    args = ["samtools", "mpileup", "-B", "-f", ref_fa, bam_file]
    logging.debug("samtools args: %s" % (' '.join(map(str,args))))
    samtools_p = subprocess.Popen(args, stdout=subprocess.PIPE)
    args = ["java", "-Xmx4g", "-jar", varscan_jar, varscan_func]
    for arg in varscan_args:
        args.extend(arg.split())
    logging.debug("varscan args: %s" % (' '.join(map(str,args))))
    outfh = open(output_file, "w")
    retcode1 = subprocess.call(args, stdin=samtools_p.stdout, stdout=outfh)
    retcode2 = samtools_p.wait()
    retcode = retcode1 + retcode2
    outfh.close()
    if retcode != 0:
        if os.path.exists(output_file):
            os.remove(output_file)
    return retcode

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("--varscan-arg", dest="varscan_args", action="append", default=[])
    parser.add_argument("varscan_jar")
    parser.add_argument("varscan_func")
    parser.add_argument("ref_fa")
    parser.add_argument("bam_file")
    parser.add_argument("output_file")
    args = parser.parse_args()
    return run_varscan(args.ref_fa, args.bam_file, args.output_file, 
                       args.varscan_jar, args.varscan_func, 
                       args.varscan_args)                     

if __name__ == '__main__': 
    sys.exit(main())