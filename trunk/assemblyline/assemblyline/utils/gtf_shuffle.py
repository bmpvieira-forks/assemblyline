'''
Created on Jan 8, 2014

@author: mkiyer
'''
'''
Created on Nov 25, 2013

@author: mkiyer
'''
import argparse
import logging
import subprocess
import collections
import os

from assemblyline.lib.base import which
from assemblyline.lib.transcript import parse_gtf
from assemblyline.lib.gtf import sort_gtf

class Interval(object):
    def __init__(self):
        self.gene_id = None
        self.chrom = None
        self.start = None
        self.end = None

def get_gene_intervals(transcripts):
    gene_map = collections.defaultdict(lambda: Interval())
    for t in transcripts:
        gene_id = t.attrs["gene_id"]
        if gene_id not in gene_map:
            g = Interval()
            g.gene_id = gene_id
            g.chrom = t.chrom
            g.start = t.start
            g.end = t.end
            gene_map[gene_id] = g
        else:
            g = gene_map[gene_id]
        # update interval
        g.start = min(g.start, t.start)
        g.end = max(g.end, t.end)
    for g in gene_map.itervalues():
        yield g
    
def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument('gtf_file')
    parser.add_argument('excl_file')
    parser.add_argument('chrom_sizes')
    parser.add_argument("output_prefix")
    args = parser.parse_args()
    prefix = args.output_prefix
    excl_file = args.excl_file
    chrom_sizes_file = args.chrom_sizes
    gtf_file = args.gtf_file
    # check command line parameters
    if which('bedtools') is None:
        parser.error('bedtools binary not found in PATH')
    if not os.path.exists(chrom_sizes_file):
        parser.error('chrom sizes file %s not found' % (chrom_sizes_file))
    gene_intervals_file = prefix + '.gene_intervals.bed'
    gene_intervals_shuffled_file = prefix + '.gene_intervals.shuffle.bed'
    shuffled_gtf_file = prefix + '.shuffle.gtf'
    sorted_shuffled_gtf_file = prefix + '.shuffle.srt.gtf'
    logging.info('Parsing GTF file')
    with open(gene_intervals_file, 'w') as f:
        for locus_transcripts in parse_gtf(open(gtf_file)):
            # find borders of locus
            locus_chrom = locus_transcripts[0].chrom
            locus_start = min(t.start for t in locus_transcripts)
            locus_end = max(t.end for t in locus_transcripts)
            logging.debug("[LOCUS] %s:%d-%d %d transcripts" % 
                          (locus_chrom, locus_start, locus_end, 
                           len(locus_transcripts)))
            for g in get_gene_intervals(locus_transcripts):
                print >>f, '\t'.join(map(str, [g.chrom, g.start, g.end, g.gene_id]))    
    # randomly shuffle genes
    logging.info("Shuffling genes")
    args = ['bedtools', 'shuffle', 
            '-excl', excl_file,
            '-i', gene_intervals_file, 
            '-g', args.chrom_sizes]
    with open(gene_intervals_shuffled_file, 'w') as fileh:
        subprocess.call(args, stdout=fileh)
    # read new gene positions
    logging.info("Reading shuffled gene intervals")
    shuffle_gene_map = {}
    with open(gene_intervals_shuffled_file) as fileh:
        for line in fileh:
            fields = line.strip().split('\t')
            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            gene_id = fields[3]
            shuffle_gene_map[gene_id] = (chrom, start, end)
    # reposition transcripts
    logging.info("Repositioning transcripts")
    with open(shuffled_gtf_file, 'w') as fileh:
        for locus_transcripts in parse_gtf(open(gtf_file)):
            # get original positions
            orig_gene_map = {}
            for g in get_gene_intervals(locus_transcripts):
                orig_gene_map[g.gene_id] = (g.chrom, g.start, g.end)
            for t in locus_transcripts:
                gene_id = t.attrs['gene_id']
                orig_chrom, orig_start, orig_end = orig_gene_map[gene_id]
                if gene_id not in shuffle_gene_map:
                    logging.warning('Gene %s [%s:%d-%d] could not be shuffled' % (gene_id, orig_chrom, orig_start, orig_end))
                    continue
                new_chrom, new_start, new_end = shuffle_gene_map[gene_id]
                # reposition transcript
                t.chrom = new_chrom
                t.start = new_start + (t.start - orig_start)
                t.end = new_start + (t.end - orig_start)
                for e in t.exons:
                    e.start = new_start + (e.start - orig_start)
                    e.end = new_start + (e.end - orig_start)
                for f in t.to_gtf_features(source='shuffle'):
                    print >>fileh, str(f)
    logging.info("Sorting GTF file")
    sort_gtf(shuffled_gtf_file, sorted_shuffled_gtf_file)

if __name__ == '__main__':
    main()