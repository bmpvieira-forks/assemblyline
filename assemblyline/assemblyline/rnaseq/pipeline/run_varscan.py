'''
Created on Feb 1, 2012

@author: oabalbin
'''
import logging
import argparse
import subprocess
import sys
import os

def run_varscan(ref_fa, 
                bam_file, 
                snv_file, 
                indel_file,
                varscan_jar,
                varscan_args):
    # run varscan snv
    args = ["samtools", "mpileup", "-f", ref_fa, bam_file]
    samtools_p = subprocess.Popen(args, stdout=subprocess.PIPE)
    args = ["java", "-jar", varscan_jar, "mpileup2snp"]
    args.extend(varscan_args)
    outfh = open(snv_file, "w")
    retcode1 = subprocess.call(args, stdin=samtools_p.stdout, stdout=outfh)
    retcode2 = samtools_p.wait()
    retcode = retcode1 + retcode2
    outfh.close()
    if retcode != 0:
        if os.path.exists(snv_file):
            os.remove(snv_file)
        return 1
    # run varscan indel
    args = ["samtools", "mpileup", "-f", ref_fa, bam_file]
    samtools_p = subprocess.Popen(args, stdout=subprocess.PIPE)    
    args = ["java", "-jar", varscan_jar, "mpileup2indel"]
    args.extend(varscan_args)
    outfh = open(indel_file, "w")
    retcode1 = subprocess.call(args, stdin=samtools_p.stdout, stdout=outfh)
    retcode2 = samtools_p.wait()
    retcode = retcode1 + retcode2    
    outfh.close()
    if retcode != 0:
        if os.path.exists(indel_file):
            os.remove(indel_file)
        return 1
    return retcode

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("--varscan-arg", dest="varscan_args", action="append", default=[])
    parser.add_argument("varscan_jar")
    parser.add_argument("ref_fa")
    parser.add_argument("bam_file")
    parser.add_argument("snv_file")
    parser.add_argument("indel_file")
    args = parser.parse_args()
    return run_varscan(args.ref_fa, args.bam_file, 
                       args.snv_file, args.indel_file,
                       args.varscan_jar, args.varscan_args)                     

if __name__ == '__main__': 
    sys.exit(main())