'''
Created on Jan 8, 2014

@author: mkiyer
'''
import argparse
import logging
import sys
import os

from assemblyline.lib.transcript import parse_gtf, NO_STRAND, POS_STRAND, NEG_STRAND, strand_int_to_str

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    # parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument('-u', type=int, dest='upstream', default=1000)
    parser.add_argument('-d', type=int, dest='downstream', default=0)
    parser.add_argument('gtf_file')
    parser.add_argument('chrom_sizes')
    args = parser.parse_args()
    upstream = args.upstream
    downstream = args.downstream
    chrom_sizes_file = args.chrom_sizes
    # check command line parameters
    if not os.path.exists(args.gtf_file):
        parser.error("GTF file %s not found" % (args.gtf_file))
    chrom_sizes = {}
    with open(chrom_sizes_file) as fileh:
        for line in fileh:
            fields = line.strip().split('\t')
            chrom_sizes[fields[0]] = int(fields[1])
    # parse
    for locus_transcripts in parse_gtf(open(args.gtf_file)):
        locus_chrom = locus_transcripts[0].chrom
        locus_start = locus_transcripts[0].start
        locus_end = max(t.end for t in locus_transcripts)
        logging.debug("[LOCUS] %s:%d-%d %d transcripts" % 
                      (locus_chrom, locus_start, locus_end, 
                       len(locus_transcripts)))
        tss_ids = set()
        for t in locus_transcripts:            
            if t.strand == NO_STRAND:
                continue
            tss_id = t.attrs['tss_id']
            if tss_id in tss_ids:
                continue
            tss_ids.add(tss_id)
            if t.strand == POS_STRAND:
                start = t.exons[0].start - upstream
                start = max(0, start)
                end = t.exons[0].start + downstream
                end = min(t.end, end)
            else:
                start = t.exons[-1].end - downstream
                start = max(t.start, start)
                end = t.exons[-1].end + upstream
                end = min(end, chrom_sizes[locus_chrom])
            print '\t'.join(map(str, [locus_chrom, start, end, tss_id, 0, strand_int_to_str(t.strand)]))

    return 0

if __name__ == "__main__":
    sys.exit(main())