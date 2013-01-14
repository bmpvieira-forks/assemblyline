'''
Created on Dec 5, 2012

@author: mkiyer
'''
import os
import sys
import logging
import argparse
import subprocess

from assemblyline.lib.gtf import GTFFeature, sort_gtf

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("reference_fasta")
    parser.add_argument("gtf_file")
    args = parser.parse_args()
    logging.info("Running cuffcompare")
    retcode = subprocess.call(["cuffcompare", "-s", args.reference_fasta, 
                               "-CG", "-r", args.gtf_file, args.gtf_file])
    if retcode != 0:
        logging.error("Error running cuffcompare")
        return retcode
    logging.info("Storing cuffcompare results")
    transcript_info_dict = {}
    for f in GTFFeature.parse(open("cuffcmp.combined.gtf")):
        p_id = f.attrs.get("p_id", None)
        tss_id = f.attrs.get("tss_id", None)
        t_id = f.attrs["oId"]
        transcript_info_dict[t_id] = (p_id, tss_id)
    logging.info("Adding attributes to original GTF")
    outfh = open("tmp.gtf", "w")
    for f in GTFFeature.parse(open(args.gtf_file)):
        t_id = f.attrs["transcript_id"]
        p_id, tss_id = transcript_info_dict[t_id]
        if p_id is not None:
            f.attrs["p_id"] = p_id
        if tss_id is not None:
            f.attrs["tss_id"] = tss_id
        print >>outfh, str(f)
    outfh.close()
    logging.info("Sorting GTF and writing 'genes.gtf'")
    sort_gtf("tmp.gtf", "genes.gtf")

if __name__ == '__main__':
    sys.exit(main())