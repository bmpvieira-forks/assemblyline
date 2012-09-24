'''
Created on Sep 6, 2012

@author: mkiyer
'''
import argparse
import logging
import os

from assemblyline.lib.gtf import GTFFeature

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("gtf_file")
    args = parser.parse_args()
    fhdict = {}
    for f in GTFFeature.parse(open(args.gtf_file)):
        sample = f.attrs["sample"]
        if not os.path.exists("%s.gtf" % sample):
            fh = open("%s.gtf" % sample, "w")
            fhdict[sample] = fh
        else:
            fh = fhdict[sample]
        print >>fh, str(f)

if __name__ == '__main__':
    main()