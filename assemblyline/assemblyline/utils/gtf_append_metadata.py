'''
Created on Jan 26, 2015

@author: mkiyer
'''
import logging
import os
import sys
import argparse

from assemblyline.lib.gtf import GTFFeature

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument('--add', dest='add', action='append', default=[])
    parser.add_argument("metadata_file")
    parser.add_argument("gtf_file")
    args = parser.parse_args()
    add_attrs = set()
    for arg in args.add:
        add_attrs.add(arg)
    add_attrs = sorted(add_attrs)

    # read metadata
    meta_dict = {}
    with open(args.metadata_file) as f:
        header_fields = f.next().strip().split('\t')
        t_id_index = header_fields.index('transcript_id')
        header_indexes = [header_fields.index(x) for x in add_attrs]
        
        for line in f:
            fields = line.strip().split('\t')
            meta_dict[fields[t_id_index]] = [fields[i] for i in header_indexes]
            #print fields[t_id_index], meta_dict[fields[t_id_index]]
    
    # read GTF
    for f in GTFFeature.parse(open(args.gtf_file)):
        t_id = f.attrs['transcript_id']
        
        if t_id in meta_dict:
            for k,v in zip(add_attrs, meta_dict[t_id]):
                f.attrs[k] = v
        else:
            for k in add_attrs:
                f.attrs[k] = 'NA'
        print str(f)
    return 0

if __name__ == '__main__':
    sys.exit(main())