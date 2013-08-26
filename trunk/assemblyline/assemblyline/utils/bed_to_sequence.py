'''
Created on Jun 22, 2012

@author: mkiyer
'''
import argparse
import logging
import sys
import os
import itertools
import string
import collections
import pysam

#Translation table for reverse Complement, with ambiguity codes
DNA_COMPLEMENT = string.maketrans( "ACGTRYKMBDHVacgtrykmbdhv", "TGCAYRMKVHDBtgcayrmkvhdb" )
RNA_COMPLEMENT = string.maketrans( "ACGURYKMBDHVacgurykmbdhv", "UGCAYRMKVHDBugcayrmkvhdb" )
#Translation table for DNA <--> RNA
DNA_TO_RNA = string.maketrans( "Tt", "Uu" )
RNA_TO_DNA = string.maketrans( "Uu", "Tt" )

def DNA_complement( sequence ):
    '''complement DNA sequence string'''
    return sequence.translate( DNA_COMPLEMENT )
def DNA_reverse_complement( sequence ):
    '''returns the reverse complement of the sequence'''
    return DNA_complement(sequence[::-1])
def to_DNA( sequence ):
    return sequence.translate( DNA_TO_RNA )
#complement RNA sequence string
def RNA_complement( sequence ):
    return sequence.translate( RNA_COMPLEMENT )
def RNA_reverse_complement( self, sequence ):
    return RNA_complement( sequence[::-1] )
def to_RNA( sequence ):
    return sequence.translate( RNA_TO_DNA )

class BEDFeature(object):    
    __slots__ = ('chrom', 'tx_start', 'tx_end', 'name', 'score', 'strand',
                 'cds_start', 'cds_end', 'exon_count', 'block_starts', 
                 'block_sizes', 'exons', 'attr_fields')

    def __str__(self):
        fields = [self.chrom,
                  str(self.tx_start),
                  str(self.tx_end),
                  self.name,
                  str(self.score),
                  self.strand,
                  str(self.cds_start),
                  str(self.cds_end),
                  '0',
                  str(self.exon_count),
                  ','.join(map(str, self.block_sizes)) + ',',
                  ','.join(map(str, self.block_starts)) + ',']
        return '\t'.join(fields)

    @staticmethod
    def from_string(line):        
        if line is None:
            return None
        line = line.strip()
        if line.startswith('#'):
            logging.debug("skipping comment line: %s" % (line))
            return None
        if line.startswith('track'):
            logging.debug("skipping track header line: %s"  % (line))
            return None
        fields = line.split('\t')
        # first six fields are required
        g = BEDFeature()
        g.chrom = fields[0]
        g.tx_start = int(fields[1])
        g.tx_end = int(fields[2])
        g.name = fields[3]
        if len(fields) <= 4:
            g.score = 0
            g.strand = '.'
        else:
            g.score = fields[4]
            g.strand = fields[5]        
        if len(fields) <= 6:
            g.cds_start = g.tx_start
            g.cds_end = g.tx_end
            g.exon_count = 1
            g.exons = [(g.tx_start, g.tx_end)]
        else:
            g.cds_start = int(fields[6])
            g.cds_end = int(fields[7])
            g.exon_count = int(fields[9])
            g.block_sizes = map(int, fields[10].split(',')[:-1])
            g.block_starts = map(int, fields[11].split(',')[:-1])            
            g.exons = []
            for start, size in itertools.izip(g.block_starts, g.block_sizes):
                g.exons.append((g.tx_start + start, g.tx_start + start + size))
        if len(fields) <= 12:
            g.attr_fields = []
        else:
            g.attr_fields = fields[12:]
        return g

    @staticmethod
    def parse(line_iter):
        for line in line_iter:
            if not line:
                continue
            if not line.strip():
                continue            
            if line.startswith("#"):
                continue
            if line.startswith("track"):
                continue
            yield BEDFeature.from_string(line)

def split_seq(seq, chars_per_line):
    pos = 0
    newseq = []
    while pos < len(seq):
        if pos + chars_per_line > len(seq):        
            endpos = len(seq)
        else:
            endpos = pos + chars_per_line
        newseq.append(seq[pos:endpos])
        pos = endpos
    return '\n'.join(newseq)

def bed_to_seq(bed_file, reference_seq_file, sep, output_fasta):
    ref_fa = pysam.Fastafile(reference_seq_file)
    for f in BEDFeature.parse(open(bed_file)):
        exon_seqs = []
        error_occurred = False
        for start, end in f.exons:
            seq = ref_fa.fetch(f.chrom, start, end)
            if (not seq) or (len(seq) < (end - start)):
                logging.warning("transcript %s exon %s:%d-%d not found in reference" % 
                                (f.name, f.chrom, start, end))
                error_occurred = True
                break
            exon_seqs.append(seq)
        if error_occurred:
            continue
        # make fasta record
        seq = sep.join(exon_seqs)
        # look for sequences containing only 'N's
        base_counts = collections.Counter(seq)
        valid_bases = sum(base_counts[x] for x in 
                          ("A","T","G","C","a","t","g","c"))
        if valid_bases == 0:
            logging.warning("transcript %s at pos %s:%d-%d lacks valid bases" %
                            (f.name, f.chrom, f.tx_start, f.tx_end))
            continue
        # reverse complement negative stranded sequences
        if f.strand == '-':
            seq = DNA_reverse_complement(seq)
        if output_fasta:
            yield ">%s %s:%d-%d[%s]" % (f.name, f.chrom, f.tx_start, f.tx_end, f.strand)
            yield seq
        else:
            yield '\t'.join([f.name, '%s:%d-%d[%s]' % (f.chrom, f.tx_start, f.tx_end, f.strand), seq])
    ref_fa.close()

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser(description="Extract transcript sequences from BED file")
    parser.add_argument("--sep", dest="sep", default='')
    parser.add_argument("--fasta", dest="fasta", action="store_true", default=False)
    parser.add_argument("ref_fasta_file", help="reference genome FASTA file")    
    parser.add_argument("bed_file")
    args = parser.parse_args()
    # check that input files exist
    if not os.path.isfile(args.ref_fasta_file):
        parser.error("Reference fasta file '%s' not found" % (args.ref_fasta_file))
    if not os.path.isfile(args.bed_file):
        parser.error("BED file '%s' not found" % (args.bed_file))
    for res in bed_to_seq(args.bed_file, args.ref_fasta_file, args.sep, args.fasta):
        print res

if __name__ == '__main__':
    sys.exit(main())
