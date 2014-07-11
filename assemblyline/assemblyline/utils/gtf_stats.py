'''
Created on Sep 5, 2012

@author: mkiyer
'''
import argparse
import logging

from assemblyline.lib.gtf import GTFFeature

def iterintrons(exons):
    e1 = exons[0]
    for j in xrange(1, len(exons)):
        e2 = exons[j]
        yield e1[1],e2[0]
        e1 = e2
        
def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("gtf_file")
    args = parser.parse_args()
    gtf_file = args.gtf_file
    genes = set()
    transcript_dict = {}
    exon_dict = {}
    for f in GTFFeature.parse(open(gtf_file)):
        if f.feature_type != 'exon':
            continue
        genes.add(f.attrs["gene_id"])
        t_id = f.attrs["transcript_id"]
        if t_id not in transcript_dict:
            transcript_dict[t_id] = (f.seqid, f.strand)
            exon_dict[t_id] = []
        exon_dict[t_id].append((f.start,f.end))

    introns = set()
    exons = set()
    for t_id in transcript_dict:
        chrom, strand = transcript_dict[t_id]
        t_exons = exon_dict[t_id]
        t_exons.sort()
        for start,end in t_exons:
            exons.add((chrom,strand,start,end))
        for start,end in iterintrons(t_exons):
            introns.add((chrom,strand,start,end))

    logging.debug("Genes: %d" % (len(genes)))
    logging.debug("Transcripts: %d" % (len(transcript_dict)))
    logging.debug("Introns: %d" % (len(introns)))
    logging.debug("Exons: %d" % (len(exons)))

if __name__ == '__main__':
    main()