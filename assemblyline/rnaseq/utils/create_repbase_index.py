'''
Created on Dec 6, 2012

@author: mkiyer
'''
import logging
import argparse
import sys
import subprocess

from assemblyline.rnaseq.lib.base import split_seq

def expand_ambiguous_bases(seq):
    """
    A     A     Adenosine     T     B
    T/U     T     Thymidine/Uridine     A     V
    G     G     Guanine     C     H
    C     C     Cytidine     G     D
    K     G or T     Keto     M     M
    M     A or C     Amino     K     K
    R     A or G     Purine     Y     Y
    Y     C or T     Pyrimidine     R     R
    S     C or G     Strong     S     W
    W     A or T     Weak     W     S
    B     C or G or T     not A (B comes after A)     V     A
    V     A or C or G     not T/U (V comes after U)     B     T/U
    H     A or C or T     not G (H comes after G)     D     G
    D     A or G or T     not C (D comes after C)     H     C
    X/N     G or A or T or C     any     N     .
    .     not G or A or T or C         .     N
    -     gap of indeterminate length             
    """
    translate = {'a': 'a',
                 't': 't',
                 'u': 't',
                 'g': 'g',
                 'c': 'c',
                 'k': 'gt',
                 'm': 'ac',
                 'r': 'ag',
                 'y': 'ct',
                 's': 'cg',
                 'w': 'at',
                 'b': 'cgt',
                 'v': 'acg',
                 'h': 'act',
                 'd': 'agt',
                 'x': 'n',
                 'n': 'n'}
    curseqs = ['']
    for c in seq:
        newseqs = []
        for curseq in curseqs: 
            for x in translate[c]:
                newseqs.append(curseq + x)
        curseqs = newseqs
    return curseqs

def create_repbase_index(input_fasta_file, index_name):
    # fix sequence names
    seqlines = []
    seqdict = {}
    for line in open(input_fasta_file):
        line = line.strip()
        if line.startswith(">"):
            if len(seqlines) > 0:
                seqdict[seqlines[0]] = seqlines[1:]
                seqlines = []
            line = line.split()[0]
        seqlines.append(line)
    output_fasta_file = index_name + ".fa"
    outfh = open(output_fasta_file, "w")
    for nameline in sorted(seqdict):
        #seqs = expand_ambiguous_bases(''.join(seqdict[nameline]))
        print >>outfh, nameline
        print >>outfh, '\n'.join(seqdict[nameline])
    outfh.close()
    # build bowtie alignment index from the fragment size exons
    logging.info("Building bowtie2 index")
    args = ["bowtie2-build", output_fasta_file, index_name]
    return subprocess.call(args)

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("input_fasta_file")
    parser.add_argument("index_name")
    args = parser.parse_args()
    return create_repbase_index(args.input_fasta_file,
                                args.index_name)

if __name__ == '__main__': 
    sys.exit(main())
