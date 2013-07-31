'''
Created on Jul 30, 2013

@author: mkiyer
'''
import logging
import argparse
import os
import sys
import collections

# project imports
import assemblyline
import assemblyline.lib.config as config
from assemblyline.lib.base import Library, GTFAttr
from assemblyline.lib.gtf import GTFFeature, sort_gtf

def read_gtf_file(filename, library_id):
    # read all transcripts
    t_dict = collections.OrderedDict()
    cur_t_id = 1
    cur_g_id = 1
    t_id_map = {}
    g_id_map = {}
    for feature in GTFFeature.parse(open(filename)):
        t_id = feature.attrs[GTFAttr.TRANSCRIPT_ID]
        # rename transcript id
        if t_id not in t_id_map:
            new_t_id = "%s.T%d" % (library_id, cur_t_id)
            t_id_map[t_id] = new_t_id
            cur_t_id += 1
        else:
            new_t_id = t_id_map[t_id]
        # rename gene id
        g_id = feature.attrs[GTFAttr.GENE_ID]
        if g_id not in g_id_map:
            new_g_id = "%s.G%d" % (library_id, cur_g_id)
            g_id_map[g_id] = new_g_id
            cur_g_id += 1
        else:
            new_g_id = g_id_map[g_id]
        # update transcript attributes
        feature.attrs[GTFAttr.TRANSCRIPT_ID] = new_t_id
        feature.attrs[GTFAttr.GENE_ID] = new_g_id
        # store feature
        if new_t_id not in t_dict:
            t_dict[new_t_id] = []
        t_dict[new_t_id].append(feature)
    return t_dict

def main():
    # setup logging
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    logging.info("AssemblyLine %s" % (assemblyline.__version__))
    logging.info("----------------------------------")
    # parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument('gtf_file_list')
    args = parser.parse_args()
    # check command line parameters
    if not os.path.exists(args.gtf_file_list):
        parser.error("gtf file list %s not found" % (args.gtf_file_list))
    # show parameters
    logging.info("Parameters:")
    logging.info("gtf file list:    %s" % (args.gtf_file_list))
    logging.info("----------------------------------")
    # parse sample table
    logging.info("Parsing GTF files")
    fileh = open(args.gtf_file_list)
    for i,line in enumerate(fileh):
        library_id = "L%d" % (i+1)
        filename = line.strip()
        t_dict = read_gtf_file(filename, library_id)
        for features in t_dict.itervalues():
            # write transcript
            for feature in features:
                print str(feature)
    logging.info("Done")
    return 0

if __name__ == '__main__':
    sys.exit(main())
