'''
Created on Aug 10, 2011

@author: mkiyer
'''
import logging
import argparse
import sys
import collections
import subprocess

# local imports
import pysam

# project imports
from assemblyline.lib.bx.intersection import Interval, IntervalTree
from assemblyline.lib.seq import DNA_reverse_complement
from assemblyline.rnaseq.lib.base import split_seq

def parse_gtf_exons(f):
    for line in f:
        fields = line.strip().split('\t')
        feature_type = fields[2]
        if feature_type != "exon":
            continue
        chrom = fields[0]
        start = int(fields[3])
        end = int(fields[4])
        strand = fields[6]
        yield chrom, start, end, strand
        
def build_exon_trees(gtf_file):
    # find all unique exons first
    chrom_exon_dict = collections.defaultdict(lambda: set())
    for chrom,start,end,strand in parse_gtf_exons(open(gtf_file)):
        chrom_exon_dict[chrom].add((start, end, strand))
    # now add exons to interval tree
    trees = collections.defaultdict(lambda: IntervalTree())
    for chrom, exons in chrom_exon_dict.iteritems():
        for start,end,strand in exons:
            trees[chrom].insert_interval(Interval(start, end, strand=strand))
    return trees

def find_unambiguous_exon_intervals(gtf_file):
    """
    returns (chrom, start, end, strand) tuples for exon
    intervals that are unique and have no overlapping
    transcripts or exons.    
    """
    trees = build_exon_trees(gtf_file) 
    for chrom,start,end,strand in parse_gtf_exons(open(gtf_file)):
        hits = set((hit.start, hit.end, hit.strand) for hit in trees[chrom].find(start, end))
        hits.add((start, end, strand))
        if len(hits) == 1:
            yield chrom, start, end, strand

def create_fragment_size_index(gtf_file, genome_fasta_file, index_name,
                               min_interval_size,
                               bowtie_build_bin):
    """
    make an alignment index containing sequences that can be used to
    assess the fragment size distribution.  these sequences must be 
    larger than the 'max_insert_size' in order to be viable for use 
    in characterizing the fragment size distribution.
    """
    # find all exons that are larger than the maximum estimated fragment size
    exons = set(coord for coord in find_unambiguous_exon_intervals(gtf_file)
                if (coord[2] - coord[1]) >= min_interval_size)
    logging.info("Found %d exons larger than %d" % (len(exons), min_interval_size))    
    # extract the nucleotide sequence of the exons
    logging.info("Extracting sequences to use for estimating the fragment "
                 "size distribution")
    ref_fa = pysam.Fastafile(genome_fasta_file)    
    frag_size_fasta_file = index_name + ".fa" 
    f = open(frag_size_fasta_file, 'w')
    for chrom, start, end, strand in exons:
        seq = ref_fa.fetch(chrom, start, end)
        if not seq:
            logging.warning("exon %s:%d-%d not found in reference" % (chrom, start, end))
            continue
        if len(seq) < (end - start):
            logging.warning("could not fetch full length seq for %s:%d-%d" % (chrom, start, end))
            continue            
        # make fasta record
        if strand == '-':
            seq = DNA_reverse_complement(seq)
        # break seq onto multiple lines
        seqlines = split_seq(seq)    
        record = (">%s:%d-%d_%s\n%s" % 
                  (chrom, start, end, strand, seqlines))
        print >>f, record
    f.close()
    ref_fa.close()
    # build bowtie alignment index from the fragment size exons
    logging.info("Building bowtie index")
    args = [bowtie_build_bin, frag_size_fasta_file, index_name]
    return subprocess.call(args)

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("--bowtie-build-bin", dest="bowtie_build_bin",
                        default="bowtie-build")
    parser.add_argument("--min-interval-size", dest="min_interval_size", 
                        type=int, default=1000)
    parser.add_argument("gtf_file")
    parser.add_argument("genome_fasta_file")
    parser.add_argument("index_name")
    args = parser.parse_args()
    return create_fragment_size_index(args.gtf_file, args.genome_fasta_file, 
                                      args.index_name,
                                      min_interval_size=args.min_interval_size,
                                      bowtie_build_bin=args.bowtie_build_bin)

if __name__ == '__main__': 
    sys.exit(main())