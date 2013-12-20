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
    parser.add_argument('--source', dest='source', default=None)
    parser.add_argument('--remove', dest='remove', action='append', default=[])
    parser.add_argument('--add', dest='add', action='append', default=[])
    parser.add_argument('gtf_file')
    args = parser.parse_args()
    gtf_file = args.gtf_file
    source = args.source
    if not os.path.exists(gtf_file):
        parser.error("GTF file '%s' not found" % (gtf_file))
    add_attrs = []
    for arg in args.add:
        k,v = arg.split(',')
        add_attrs.append((k,v))
    rm_attrs = []
    for arg in args.remove:
        rm_attrs.append(arg)
    i = 0
    for f in GTFFeature.parse(open(args.gtf_file)):
        if source is not None:
            f.source = source
        for k,v in add_attrs:
            if k in f.attrs:
                parser.error('attribute %s already in feature' % (k))
            f.attrs[k] = v
        for k in rm_attrs:
            if k in f.attrs:
                del f.attrs[k]
        print str(f)
        i += 1
        if (i % 100000) == 0:
            logging.debug('finished %d lines' % (i))
    return 0

if __name__ == '__main__':
    sys.exit(main())
