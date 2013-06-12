'''
Created on Jan 19, 2013

@author: mkiyer
'''
import logging
import argparse
import collections
import sys
import os

from assemblyline.lib.bx.intersection import Interval, IntervalTree
from assemblyline.lib.bed import BEDFeature
from assemblyline.lib.gtf import GTFFeature, parse_loci

def build_interval_tree_from_bed(bed_file):
    trees = collections.defaultdict(lambda: IntervalTree())
    for f in BEDFeature.parse(open(bed_file)):
        tree = trees[f.chrom]
        for start,end in f.exons:
            tree.insert_interval(Interval(start, end, strand=f.strand, value=f.name))
    return trees

def annotate_gtf(gtf_file, bed_dbs):
    # read reference databases
    bed_trees = []
    for name,filename in bed_dbs:
        logging.debug("Loading BED db '%s' file '%s'" % (name,filename))
        trees = build_interval_tree_from_bed(filename)
        bed_trees.append((name, trees))
    # parse gtf file and annotate
    logging.debug("Annotating GTF")
    for lines in parse_loci(open(gtf_file)):
        features = []
        transcripts = []
        transcript_matches = collections.defaultdict(lambda: collections.defaultdict(lambda: set()))
        for line in lines:
            f = GTFFeature.from_string(line)
            features.append(f)
            t_id = f.attrs['transcript_id']
            if f.feature_type == 'transcript':
                transcripts.append(f)
            elif f.feature_type == 'exon':
                for dbname,dbtrees in bed_trees:
                    # intersect this exon with features
                    hits = dbtrees[f.seqid].find(f.start, f.end)                        
                    matches = set(hit.value for hit in hits if hit.strand == f.strand)
                    f.attrs[dbname] = ','.join(sorted(matches))
                    # update transcript level matches
                    transcript_matches[t_id][dbname].update(matches)
        # set transcript annotations
        for f in transcripts:
            t_id = f.attrs['transcript_id']
            for dbname,dbtrees in bed_trees:
                matches = transcript_matches[t_id][dbname]
                f.attrs[dbname] = ','.join(sorted(matches))
        # write features
        for f in features:
            print str(f)
    logging.debug("Done")

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("--bed", dest="bed_args", action="append", default=[])
    parser.add_argument("gtf_file")
    args = parser.parse_args()
    if not os.path.exists(args.gtf_file):
        parser.error("GTF file %s not found" % (args.gtf_file))
    bed_dbs = []
    for arg in args.bed_args:
        name,filename = arg.split(",")
        bed_dbs.append((name,filename))
    annotate_gtf(args.gtf_file, bed_dbs)
    return 0

if __name__ == '__main__':
    sys.exit(main())