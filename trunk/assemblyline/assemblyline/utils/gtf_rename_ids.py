'''
Created on Jul 31, 2013

@author: mkiyer
'''
import logging
import argparse
import os
import sys

# project imports
from assemblyline.lib.gtf import GTFFeature

def main():
    # setup logging
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    # parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument('gtf_file')
    args = parser.parse_args()
    # check command line parameters
    if not os.path.exists(args.gtf_file):
        parser.error("gtf file %s not found" % (args.gtf_file))
    cur_t_id = 1
    cur_g_id = 1
    cur_tss_id = 1
    t_id_map = {}
    g_id_map = {}
    tss_id_map = {}
    for feature in GTFFeature.parse(open(args.gtf_file)):
        t_id = feature.attrs['transcript_id']
        g_id = feature.attrs['gene_id']
        tss_id = feature.attrs['tss_id']
        if t_id not in t_id_map:
            new_t_id = "T%06d" % (cur_t_id)
            t_id_map[t_id] = new_t_id
            cur_t_id += 1
        else:
            new_t_id = t_id_map[t_id]
        if g_id not in g_id_map:
            new_g_id = "G%06d" % (cur_g_id)
            g_id_map[g_id] = new_g_id
            cur_g_id += 1
        else:
            new_g_id = g_id_map[g_id]
        if tss_id not in tss_id_map:
            new_tss_id = "TSS%06d" % (cur_tss_id)
            tss_id_map[tss_id] = new_tss_id
            cur_tss_id += 1
        else:
            new_tss_id = tss_id_map[tss_id]
        # update transcript attributes
        new_attrs = {'transcript_id': new_t_id,
                     'gene_id': new_g_id,
                     'tss_id': new_tss_id}
        if 'exon_number' in feature.attrs:
            new_attrs['exon_number'] = feature.attrs['exon_number']
        feature.attrs = new_attrs
        print feature
    return 0

if __name__ == '__main__':
    sys.exit(main())
