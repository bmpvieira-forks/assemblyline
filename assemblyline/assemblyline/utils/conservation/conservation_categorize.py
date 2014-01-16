'''
Created on Jan 7, 2014

@author: mkiyer
'''
import logging
import argparse
import os
import sys
import collections

import assemblyline
from assemblyline.lib.gtf import GTFFeature, sort_gtf
from assemblyline.lib.transcript import parse_gtf, strand_int_to_str, NO_STRAND, POS_STRAND, NEG_STRAND
from assemblyline.lib.assemble.transcript_graph import \
    find_exon_boundaries, split_exons

def categorize(transcripts):
    # find the intron domains of the transcripts
    boundaries = find_exon_boundaries(transcripts)
    cds_sense = set()
    cds_unstranded = set()
    test_transcripts = []
    for t in transcripts:
        if 'cds' in t.attrs:            
            for n in split_exons(t, boundaries):
                cds_sense.add((t.strand, n[0], n[1]))
                cds_unstranded.add(n)
        else:
            test_transcripts.append(t)
    node_dict = {}
    noncoding_transcripts = []
    utr_unstranded = set()
    for t in test_transcripts:
        mcat = t.attrs['category']
        tcat = t.attrs['transcript_category']
        has_cds = False
        nodes = []
        for n in split_exons(t, boundaries):
            nstrand = (t.strand, n[0], n[1])
            nodes.append(nstrand)
            if ((t.strand == NO_STRAND) and (n in cds_unstranded)):
                node_dict[nstrand] = (tcat, mcat, 'cds')
                has_cds = True
            elif nstrand in cds_sense:
                node_dict[nstrand] = (tcat, mcat, 'cds')
                has_cds = True
        if has_cds:
            for n in nodes:
                if n not in node_dict:
                    node_dict[nstrand] = (tcat, mcat, 'utr')
                    utr_unstranded.add((n[1],n[2]))
        else:
            noncoding_transcripts.append(t)

    for t in noncoding_transcripts:
        mcat = t.attrs['category']
        tcat = t.attrs['transcript_category']
        for n in split_exons(t, boundaries):
            nstrand = (t.strand, n[0], n[1])
            if nstrand not in node_dict:
                # check antisense
                if n in cds_unstranded:
                    node_dict[nstrand] = (tcat, mcat, 'cds_antisense')
                elif n in utr_unstranded:
                    node_dict[nstrand] = (tcat, mcat, 'utr_antisense')
                else:
                    node_dict[nstrand] = (tcat, mcat, 'noncoding')
    for n,d in node_dict.iteritems():
        strand, start, end = n
        tcat, mcat, ccat = d
        yield start, end, strand, tcat, mcat, ccat

def get_cds_features(gtf_file):
    cds = collections.defaultdict(lambda: set())
    i = 0
    for f in GTFFeature.parse(open(gtf_file)):
        if f.feature_type == "CDS":
            cds[f.seqid].add((f.start, f.end, f.strand))
        i += 1
        if (i % 100000) == 0:
            logging.debug("Parsed %d features" % (i))
    logging.debug("Returning CDS transcripts")
    t_id = 1
    for chrom in sorted(cds):
        for start,end,strand in sorted(cds[chrom]):
            for feature_type in ('transcript', 'exon'):
                f = GTFFeature()
                f.seqid = chrom
                f.source = 'cds'
                f.feature_type = feature_type
                f.start = start
                f.end = end
                f.score = 0
                f.strand = strand
                f.phase = '.'
                f.attrs = {'cds': 1, 
                           'transcript_id': 'CDS%08d' % (t_id)}
                yield f
            t_id += 1

def main():
    # parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("-v", "--verbose", action="store_true", 
                        dest="verbose", default=False)
    parser.add_argument("ref_gtf_file")
    parser.add_argument("gtf_file")
    args = parser.parse_args()
    # set logging level
    if args.verbose:
        level = logging.DEBUG
    else:
        level = logging.INFO
    logging.basicConfig(level=level,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    # check command line parameters
    if not os.path.exists(args.ref_gtf_file):
        parser.error("GTF file %s not found" % (args.ref_gtf_file))
    if not os.path.exists(args.gtf_file):
        parser.error("GTF file %s not found" % (args.gtf_file))
    logging.info("AssemblyLine %s" % (assemblyline.__version__))
    logging.info("----------------------------------")   
    # show parameters
    logging.info("Parameters:")
    logging.info("verbose logging:       %s" % (args.verbose))
    logging.info("ref gtf file:          %s" % (args.ref_gtf_file))
    logging.info("assembly gtf file:     %s" % (args.gtf_file))
    # find CDS regions
    if not os.path.exists('tmp.srt.gtf'):
        with open('tmp.gtf', 'w') as outfileh:
            logging.info("Reading CDS regions from reference GTF")
            for f in get_cds_features(args.ref_gtf_file):
                print >>outfileh, str(f)
            logging.info("Reading transcripts from assembly GTF")
            i = 0
            for f in GTFFeature.parse(open(args.gtf_file)):
                print >>outfileh, str(f)
                i += 1
                if i % 100000 == 0:
                    logging.debug("Parsed %d transcripts" % (i))
        logging.info("Sorting GTF file")
        sort_gtf('tmp.gtf', 'tmp.srt.gtf')
    for locus_transcripts in parse_gtf(open('tmp.srt.gtf')):
        locus_chrom = locus_transcripts[0].chrom
        locus_start = locus_transcripts[0].start
        locus_end = max(t.end for t in locus_transcripts)
        logging.debug("[LOCUS] %s:%d-%d %d transcripts" % 
                      (locus_chrom, locus_start, locus_end, 
                       len(locus_transcripts)))
        for start, end, strand, m, t, c in categorize(locus_transcripts):
            fields = [locus_chrom, str(start), str(end), '%s|%s|%s' % (m,t,c), '0', strand_int_to_str(strand)]
            print '\t'.join(fields)    
    return 0

if __name__ == "__main__":
    sys.exit(main())