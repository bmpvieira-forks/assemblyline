'''
Created on Aug 17, 2012

@author: mkiyer
'''
import logging
import argparse

from assemblyline.lib.transcript_parser import parse_gtf

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("--length", type=int, dest="length", default=200)
    parser.add_argument("--utr-length", type=int, dest="utr_length", default=25)
    parser.add_argument("--antisense-exons", type=int, dest="antisense_exons", default=2)
    parser.add_argument("--intronic-exons", type=int, dest="intronic_exons", default=2)
    parser.add_argument("--intergenic-exons", type=int, dest="intergenic_exons", default=2)
    parser.add_argument("--intergenic-dist", type=int, dest="intergenic_dist", default=1000)
    parser.add_argument("gtf_file")
    args = parser.parse_args()    
    lncrna_categories = set(("intronic", "intergenic", "antisense", "ncrna"))
    logging.debug("Retrieving lncrna features from GTF")
    total_transcripts = 0
    lncrnas = 0
    intergenic = 0
    intronic = 0
    antisense = 0
    for transcripts in parse_gtf(open(args.gtf_file)):
        total_transcripts += len(transcripts)
        for t in transcripts:
            category = t.attrs["category"]
            if category not in lncrna_categories:
                continue
            # check length requirement
            if t.length < args.length:
                continue
            if category == "intergenic":
                dist = int(t.attrs["nearest_dist"])
                # exclude "intergenic" lncrnas that are close to known genes
                if dist < args.intergenic_dist:
                    continue
                if len(t.exons) < args.intergenic_exons:
                    continue
                intergenic += 1
            elif category == "intronic":
                if len(t.exons) < args.intronic_exons:
                    continue
                intronic += 1
            elif category == "antisense":
                if len(t.exons) < args.antisense_exons:
                    continue
                antisense += 1
            # check utr exon length requirement
            #if len(t.exons) == 2:
            #    if (t.exons[0].end - t.exons[0].start) < args.utr_length:
            #        continue
            #    if (t.exons[-1].end - t.exons[-1].start) < args.utr_length:
            #        continue
            # output
            for f in t.to_gtf_features():
                print str(f)
            lncrnas += 1
    logging.debug("Read %d lncrna transcripts out of %d total transcripts" % (lncrnas, total_transcripts))
    logging.debug("intergenic: %d" % (intergenic))
    logging.debug("intronic: %d" % (intronic))
    logging.debug("antisense: %d" % (antisense))

if __name__ == '__main__':
    main()