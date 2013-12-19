'''
Created on Dec 18, 2013

@author: mkiyer
'''
import os
import sys
import logging
import argparse

from assemblyline.lib.gtf import GTFFeature

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument('gtf_file')
    parser.add_argument('key')
    parser.add_argument('value')
    args = parser.parse_args()
    gtf_file = args.gtf_file
    key = args.key
    value = args.value
    if not os.path.exists(gtf_file):
        parser.error("GTF file '%s' not found" % (gtf_file))
    for f in GTFFeature.parse(open(args.gtf_file)):
        if key in f.attrs:
            parser.error('attribute %s already in feature' % (key))
        f.attrs[key] = value
        print str(f)
    return 0

if __name__ == '__main__':
    sys.exit(main())
