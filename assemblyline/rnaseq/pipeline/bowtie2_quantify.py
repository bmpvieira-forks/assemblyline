'''
Created on Dec 6, 2012

@author: mkiyer
'''
import argparse
import logging
import sys
import os
import array

# local import
import pysam
# project imports
import assemblyline.rnaseq.lib.config as config

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("sam_file")
    parser.add_argument("output_file")
    args = parser.parse_args()
    # check args
    if not os.path.exists(args.sam_file):
        parser.error("sam file %s not found" % (args.sam_file))
    # open files
    suffix = os.path.splitext(args.sam_file)[-1]
    mode = 'rb' if suffix == '.bam' else 'r'
    samfh = pysam.Samfile(args.sam_file, mode)
    nrefs = len(samfh.references)
    # init counts
    def make_array(n):
        return array.array('L', (0 for x in xrange(n)))
    # Value of UU indicates the read was not part of a pair. 
    # Value of CP indicates the read was part of a pair and the pair aligned concordantly. 
    # Value of DP indicates the read was part of a pair and the pair aligned discordantly. 
    # Value of UP indicates the read was part of a pair but the pair failed to aligned either concordantly or discordantly.
    counts = {'UU': make_array(nrefs),
              'CP': make_array(nrefs),
              'DP': make_array(nrefs),
              'UP': make_array(nrefs)}
    unmapped = 0
    for r in samfh:
        if r.is_unmapped:
            unmapped += 1
        else:
            altype = r.opt('YT')
            counts[altype][r.tid] += 1
    # print output
    outfh = open(args.output_file, "w")
    print >>outfh, '\t'.join(["#rname", "length", "paired", 
                              "paired_discordant", "paired_orphan", 
                              "unpaired", "total"])
    for i in xrange(len(samfh.references)):
        paired = counts['CP'][i]
        paired_discordant = counts['DP'][i]
        paired_orphan = counts['UP'][i]
        unpaired = counts['UU'][i]
        total = paired + paired_discordant + paired_orphan + unpaired
        fields = [samfh.references[i],
                  samfh.lengths[i],
                  paired,
                  paired_discordant,
                  paired_orphan,
                  unpaired,
                  total]
        print >>outfh, '\t'.join(map(str, fields))
    logging.info("found %d unmapped reads" % (unmapped))
    return config.JOB_SUCCESS

if __name__ == '__main__': 
    sys.exit(main())