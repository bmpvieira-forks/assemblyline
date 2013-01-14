'''
Created on Sep 5, 2012

@author: mkiyer
'''
import argparse
import logging

from assemblyline.lib.gtf import GTFFeature

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("gtf_file")
    args = parser.parse_args()
    genes = set()
    transcripts = set()
    exons = set()
    for f in GTFFeature.parse(open(args.gtf_file)):
        genes.add(f.attrs["gene_id"])
        transcripts.add(f.attrs["transcript_id"])
        exons.add((f.strand,f.start,f.end))
    logging.debug("Genes: %d" % (len(genes)))
    logging.debug("Transcripts: %d" % (len(transcripts)))
    logging.debug("Exons: %d" % (len(exons)))

if __name__ == '__main__':
    main()