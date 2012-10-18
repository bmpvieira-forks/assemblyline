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
    parser.add_argument("--gtf-split-attr", dest="gtf_split_attr", default="library_id")
    parser.add_argument("gtf_file")
    args = parser.parse_args()
    fhdict = {}
    gtf_split_attr = args.gtf_split_attr
    for f in GTFFeature.parse(open(args.gtf_file)):
        if gtf_split_attr not in f.attrs:
            val = "na_missing"
        else:
            val = f.attrs[args.gtf_split_attr]
        filename = "%s.gtf" % (val)
        if not os.path.exists(filename):
            fh = open(filename, "w")
            fhdict[val] = fh
        else:
            fh = fhdict[val]
        print >>fh, str(f)
    for fh in fhdict.itervalues():
        fh.close()

if __name__ == '__main__':
    main()