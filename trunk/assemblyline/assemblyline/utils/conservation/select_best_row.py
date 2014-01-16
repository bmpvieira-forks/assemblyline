'''
Created on Jan 4, 2014

@author: mkiyer
'''
import sys
import argparse
import logging
import collections

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument('--groupby', dest='groupby', default='gene_id')
    parser.add_argument('-f', '--field', dest='field', default='mean')
    parser.add_argument('metadata_file')
    args = parser.parse_args()
    metadata_file = args.metadata_file
    value_field = args.field
    groupby_field = args.groupby
    with open(metadata_file) as f:
        header_fields = f.next().strip().split('\t')
        print '\t'.join(header_fields)
        group_ind = header_fields.index(groupby_field)
        value_ind = header_fields.index(value_field)
        d = collections.defaultdict(lambda: [])
        for line in f:
            fields = line.strip().split('\t')
            d[fields[group_ind]].append(line)
        for g_id in sorted(d):
            lines = d[g_id]
            best_value = None
            best_fields = None
            for line in lines:
                fields = line.strip().split('\t')
                v = fields[value_ind]
                if v == 'None' or v == 'NA' or v == 'nan':
                    v = None
                else:
                    v = float(v)
                if (best_value is None) or (v is None) or (best_value < v):
                    best_value = v
                    best_fields = fields
            print '\t'.join(best_fields)   

if __name__ == '__main__':
    sys.exit(main())
