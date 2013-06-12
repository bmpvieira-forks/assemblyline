'''
Created on Aug 17, 2012

@author: mkiyer
'''
import logging
import argparse

from assemblyline.lib.transcript import parse_gtf

IGNORE_SOURCES = set(["IG_C_gene",
                      "IG_C_pseudogene",
                      "IG_D_gene",
                      "IG_J_gene",
                      "IG_J_pseudogene",
                      "IG_V_gene",
                      "IG_V_pseudogene",
                      "Mt_rRNA",
                      "Mt_tRNA",
                      "rRNA",
                      "TR_C_gene",
                      "TR_D_gene",
                      "TR_J_gene",
                      "TR_J_pseudogene",
                      "TR_V_gene",
                      "TR_V_pseudogene"])

PSEUDOGENE_SOURCES = set(["processed_pseudogene",
                          "processed_transcript",
                          "polymorphic_pseudogene",
                          "pseudogene",
                          "retrotransposed",
                          "transcribed_processed_pseudogene",
                          "transcribed_unprocessed_pseudogene",
                          "unitary_pseudogene",
                          "unprocessed_pseudogene"])

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("--length", type=int, dest="length", default=250)
    parser.add_argument("--ncrna-exons", type=int, dest="ncrna_exons", default=1)
    parser.add_argument("--pseudogene-exons", type=int, dest="pseudogene_exons", default=1)
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
    pseudogene = 0
    ncrna = 0
    for transcripts in parse_gtf(open(args.gtf_file)):
        total_transcripts += len(transcripts)
        for t in transcripts:
            # throw out protein coding genes
            category = t.attrs["category"]
            if category not in lncrna_categories:
                continue
            # throw out transcripts that overlap certain classes of 
            # transcripts
            annotation_sources = set(t.attrs["annotation_sources"].split(','))
            if not annotation_sources.isdisjoint(IGNORE_SOURCES):
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
            else:
                # keep multi-exonic pseudogenes
                if not annotation_sources.isdisjoint(PSEUDOGENE_SOURCES):
                    if len(t.exons) < args.pseudogene_exons:
                        continue
                    pseudogene += 1
                else:
                    if len(t.exons) < args.ncrna_exons:
                        continue
                    ncrna += 1
            # output
            for f in t.to_gtf_features():
                print str(f)
            lncrnas += 1
    logging.debug("Read %d lncrna transcripts out of %d total transcripts" % (lncrnas, total_transcripts))
    logging.debug("intergenic: %d" % (intergenic))
    logging.debug("intronic: %d" % (intronic))
    logging.debug("antisense: %d" % (antisense))
    logging.debug("pseudogene: %d" % (pseudogene))
    logging.debug("ncrna: %d" % (ncrna))

if __name__ == '__main__':
    main()