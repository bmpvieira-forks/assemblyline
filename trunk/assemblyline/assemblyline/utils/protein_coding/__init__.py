'''
Created on Jan 8, 2014

@author: mkiyer
'''
import argparse
import logging
import subprocess
import collections
import os
import pysam

from assemblyline.lib.base import which
from assemblyline.lib.transcript import parse_gtf

PFAM_FILES = ['Pfam-A.hmm',
              'Pfam-A.hmm.dat',
              'Pfam-A.hmm.h3f',
              'Pfam-A.hmm.h3i',
              'Pfam-A.hmm.h3m',
              'Pfam-A.hmm.h3p',
              'Pfam-B.hmm',
              'Pfam-B.hmm.dat',
              'Pfam-B.hmm.h3f',
              'Pfam-B.hmm.h3i',
              'Pfam-B.hmm.h3m',
              'Pfam-B.hmm.h3p']
def check_pfam_dir(path):
    for pfam_file in PFAM_FILES:
        if not os.path.exists(os.path.join(path, pfam_file)):
            return False
    return True
    
def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument('--pfam', dest='pfam_dir')
    parser.add_argument('--genome-fasta', dest='genome_fasta_file')
    parser.add_argument('gtf_file')
    args = parser.parse_args()
    # get args
    pfam_dir = args.pfam_dir
    genome_fasta_file = args.genome_fasta_file
    gtf_file = args.gtf_file
    # check command line parameters
    if which('pfamscan.pl') is None:
        parser.error("'pfamscan.pl' not found in PATH")
    if which('hmmscan') is None:
        parser.error("'hmmscan' not found in PATH")
    if which('signalp') is None:
        parser.error("'signalp' not found in PATH")
    if not check_pfam_dir(pfam_dir):
        parser.error("Required Pfam-A and Pfam-B files not found at '%s'" % (pfam_dir))
    if not os.path.exists(genome_fasta_file):
        parser.error("Genome FASTA file '%s' not found" % (genome_fasta_file))
    if not os.path.exists(gtf_file):
        parser.error("GTF file '%s' not found" % (gtf_file))

    return

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