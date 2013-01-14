'''
Created on Dec 19, 2012

@author: mkiyer
'''
import logging
import argparse
import sys
import os
import subprocess

def bed_to_bigbed(input_file, output_file, 
                  chrom_sizes_file,
                  score_to_name, tmp_dir):
    prefix = os.path.splitext(os.path.basename(input_file))[0]
    raw_bed_file = os.path.join(tmp_dir, prefix + ".bed")
    sorted_bed_file = os.path.join(tmp_dir, prefix + ".sorted.bed")
    # remove track lines and fix score
    logging.debug("Removing track lines and fixing score field")
    try:
        fileh = open(raw_bed_file, "w")
        for line in open(input_file):
            if line.startswith("#"):
                continue
            if line.startswith("track"):
                continue
            fields = line.strip().split('\t')                
            if len(fields) >= 4:
                if fields[4] == ".":
                    fields[4] = "0"
                else:
                    if int(fields[4]) > 1000:
                        fields[4] = "1000"
                    elif int(fields[4]) < 0:
                        fields[4] = "0"
                if score_to_name:
                    fields[3] = "%s[%s]" % (fields[3],fields[4]) 
            print >>fileh, '\t'.join(fields)
        fileh.close()
    except Exception:
        if os.path.exists(raw_bed_file):
            os.remove(raw_bed_file)
        return 1
    # sort file
    logging.debug("Sorting bed file")
    args = ['sort', '-k1,1', '-k2,2n', raw_bed_file]
    fileh = open(sorted_bed_file, "w")
    retcode = subprocess.call(args, stdout=fileh)
    fileh.close()
    if retcode != 0:
        if os.path.exists(raw_bed_file):
            os.remove(raw_bed_file)
        if os.path.exists(sorted_bed_file):
            os.remove(sorted_bed_file)
        return 1
    # run the bedToBigBed program 
    logging.debug("Making bigBed file")
    args = ["bedToBigBed", sorted_bed_file, chrom_sizes_file, output_file]
    retcode = subprocess.call(args)
    if os.path.exists(raw_bed_file):
        os.remove(raw_bed_file)
    if os.path.exists(sorted_bed_file):
        os.remove(sorted_bed_file)
    if retcode != 0:
        if os.path.exists(output_file):
            os.remove(output_file)
        return 1
    return 0

def main():
    # setup logging
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    # command line parsing
    parser = argparse.ArgumentParser()
    parser.add_argument('--score-to-name', dest="score_to_name", 
                        action="store_true", default=False)
    parser.add_argument('--tmp-dir', dest="tmp_dir", default="/tmp")
    parser.add_argument('bed_file')
    parser.add_argument('chrom_sizes_file')
    parser.add_argument('bigbed_file')
    args = parser.parse_args()
    return bed_to_bigbed(args.bed_file, args.bigbed_file, 
                         args.chrom_sizes_file,
                         score_to_name=args.score_to_name, 
                         tmp_dir=args.tmp_dir)

if __name__ == '__main__':
    sys.exit(main())