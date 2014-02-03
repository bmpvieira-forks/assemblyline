'''
Created on Feb 2, 2014

@author: mkiyer
'''
import logging
import argparse
import os
import sys

import pysam

from assemblyline.lib.transcript import parse_gtf, NEG_STRAND, NO_STRAND, POS_STRAND, strand_int_to_str
from protein_coding_potential import get_transcript_dna_sequence, GENETIC_CODE, transcript_to_genome_pos, genome_interval_to_exons, ORFInfo

def orf_to_genome(t, start, end):
    reverse = (t.strand == NEG_STRAND)
    g_orf_start = transcript_to_genome_pos(t, start, reverse)
    g_orf_end = transcript_to_genome_pos(t, end, reverse)
    if reverse:
        g_orf_start, g_orf_end = g_orf_end, g_orf_start    
    g_orf_end += 1 # end position is inclusive so adjust
    orf_exons = genome_interval_to_exons(g_orf_start, g_orf_end, t)
    return g_orf_start, g_orf_end, orf_exons

def translate_orf(seq):
    aa_seq = []
    for n in range(0, len(seq), 3):
        codon = seq[n:n+3]
        aa = GENETIC_CODE.get(codon, '.')
        aa_seq.append(aa)
        if aa == '*':
            break
    aa_seq = ''.join(aa_seq).rstrip('.')
    return aa_seq

def find_first_orf(t, ref_fa):
    orf = ORFInfo()
    orf.transcript_id = t.attrs['transcript_id']
    orf.gene_id = t.attrs['gene_id']
    orf.chrom = t.chrom
    # get transcript sequence
    seq = get_transcript_dna_sequence(t, ref_fa)
    # find first ATG in sequence
    start = seq.find('ATG')
    if start == -1:
        orf.start = t.start
        orf.end = t.start
        orf.strand = '.'
        orf.exons = []
        orf.seq = ''
    else:
        aa_seq = translate_orf(seq[start:])
        end = start + 3 * len(aa_seq) - 1
        orf_start, orf_end, orf_exons = \
            orf_to_genome(t, start, end)
        orf.start = orf_start
        orf.end = orf_end
        orf.strand = strand_int_to_str(t.strand)
        orf.exons = orf_exons
        orf.seq = aa_seq
    return orf

def find_first_orfs(gtf_file, genome_fasta_file):
    #
    # extract transcript DNA sequences, translate to protein, and
    # search for first ORF
    #
    logging.debug('Finding ORFs in transcript sequences')
    # open genome fasta file
    ref_fa = pysam.Fastafile(genome_fasta_file)
    num_finished = 1
    for locus_transcripts in parse_gtf(open(gtf_file)):
        for t in locus_transcripts:
            if t.strand == NO_STRAND:
                t.strand = POS_STRAND
                orfpos = find_first_orf(t, ref_fa)
                t.strand = NEG_STRAND
                orfneg = find_first_orf(t, ref_fa)
                if len(orfpos.seq) >= len(orfneg.seq):
                    yield orfpos
                else:
                    yield orfneg 
            else:
                yield find_first_orf(t, ref_fa)
            if (num_finished % 10000) == 0:
                logging.debug('Processed %d transcripts' % (num_finished))
            num_finished += 1
    # cleanup
    ref_fa.close()
    return

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument('--genome-fasta', dest='genome_fasta_file')
    parser.add_argument('gtf_file')
    parser.add_argument('output_prefix')
    args = parser.parse_args()
    # get args
    genome_fasta_file = args.genome_fasta_file
    gtf_file = args.gtf_file
    output_prefix = args.output_prefix
    # check command line parameters
    if not os.path.exists(genome_fasta_file):
        parser.error("Genome FASTA file '%s' not found" % (genome_fasta_file))
    if not os.path.exists(gtf_file):
        parser.error("GTF file '%s' not found" % (gtf_file))
    bedfh = open(output_prefix + '.bed', 'w')
    outfh = open(output_prefix + '.txt', 'w')
    for orf in find_first_orfs(gtf_file, genome_fasta_file):
        if len(orf.seq) == 0:
            continue
        else:
            print >>bedfh, '\t'.join(orf.to_bed(orf.transcript_id))
        print >>outfh, '\t'.join([orf.transcript_id, str(len(orf.seq))])
    bedfh.close()
    outfh.close()
    return 0


if __name__ == '__main__':
    sys.exit(main())