'''
Created on Dec 10, 2012

@author: mkiyer
'''
import logging
import argparse
import sys
import os
import subprocess

import assemblyline.rnaseq.lib.picard as picard
import assemblyline.rnaseq.lib.config as config
from assemblyline.rnaseq.lib.libtable import FR_UNSTRANDED
from assemblyline.rnaseq.lib.inspect import RnaseqLibraryCharacteristics
from assemblyline.rnaseq.lib.base import up_to_date, many_up_to_date

STRAND_COLOR_MAP = {".": "0,128,0",
                    "+": "255,0,0",
                    "-": "0,0,255"}

def make_coverage_map(bam_file, bigwig_file, strand, scale_factor, 
                      chrom_sizes_file, tmp_dir):
    # generate bedgraph coverage file
    logging.info("Generating unstranded coverage bigwig file")
    args = ["bedtools", "genomecov", "-bg", "-split", 
            "-scale", scale_factor,
            "-ibam", bam_file]
    if strand != ".":
        args.extend(["-strand", strand])
    logging.debug("bedtools args: %s" % (map(str, args)))
    bedgraph_file = os.path.join(tmp_dir, "coverage.bedgraph")
    outfh = open(bedgraph_file, "w")
    retcode = subprocess.call(map(str, args), stdout=outfh)
    outfh.close()
    if retcode != 0:
        logging.error("Error during bedtools genomecov")
        if os.path.exists(bedgraph_file):
            os.remove(bedgraph_file)
        return 1
    # convert bedgraph to bigwig coverage file
    logging.info("Create bigWig file to display coverage")        
    args = ["bedGraphToBigWig", bedgraph_file, chrom_sizes_file, bigwig_file]  
    logging.debug("bedGraphToBigWig args: %s" % (map(str, args)))
    retcode = subprocess.call(args)
    if retcode != 0:
        logging.error("Error during bedGraphToBigWig")
        if os.path.exists(bigwig_file):
            os.remove(bigwig_file)
        return 1
    if os.path.exists(bedgraph_file):
        os.remove(bedgraph_file)
    return 0

def main():
    # setup logging
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    # command line parsing
    parser = argparse.ArgumentParser()
    parser.add_argument('--big-data-url', dest="big_data_url", default="")
    parser.add_argument('--track-name', dest="track_name", default="coverage")
    parser.add_argument('--track-desc', dest="track_desc", default="coverage")
    parser.add_argument('--tmp-dir', dest="tmp_dir", default="/tmp")
    parser.add_argument('--scale', dest='scale', action="store_true", default=False)
    parser.add_argument('library_metrics_file')
    parser.add_argument('picard_alignment_metrics_file')
    parser.add_argument('chrom_sizes_file')
    parser.add_argument('bam_file')
    parser.add_argument('bigwig_prefix')
    parser.add_argument('ucsc_track_file')
    args = parser.parse_args()
    # determine normalization factor from picard metrics
    scale_factor = 1.0 
    if args.scale:
        logging.debug("Reading library size from picard alignment metrics file")
        library_size = picard.get_total_reads(args.picard_alignment_metrics_file)
        logging.debug("Found library size: %d" % (library_size))
        if library_size > 0:
            scale_factor = 1.0e6 / float(library_size)
    # predict library type
    logging.debug("Predicting library type from metrics file")
    obj = RnaseqLibraryCharacteristics.from_file(open(args.library_metrics_file))
    predicted_library_type = obj.predict_library_type(config.STRAND_SPECIFIC_CUTOFF_FRAC)
    logging.debug("Predicted library type: %s" % (predicted_library_type))
    if predicted_library_type == FR_UNSTRANDED:
        strand_files = [(".", args.bigwig_prefix + ".bw")]
    else:
        strand_files = [("+", args.bigwig_prefix + "_pos.bw"),
                        ("-", args.bigwig_prefix + "_neg.bw")] 
    for strand,bigwig_file in strand_files:
        msg = "Generating coverage map for strand=%s bigwig_file=%s" % (strand,bigwig_file)
        if up_to_date(bigwig_file, args.bam_file):
            logging.info("[SKIPPED] %s" % (msg))
        else:
            logging.info(msg)
            retcode = make_coverage_map(args.bam_file, bigwig_file, strand, 
                                        scale_factor, args.chrom_sizes_file, 
                                        args.tmp_dir)
            if retcode != 0:
                logging.error("Error")
                if os.path.exists(bigwig_file):
                    os.remove(bigwig_file)
                return 1
    # make ucsc track file
    msg = "Generating UCSC track files"
    if up_to_date(args.ucsc_track_file, args.bam_file):
        logging.info("[SKIPPED] %s" % (msg))
    else:
        logging.info(msg)
        fileh = open(args.ucsc_track_file, "w")
        for strand,bigwig_file in strand_files:
            ucsc_url = "%s%s" % (args.big_data_url, bigwig_file)
            track_line = ' '.join(['track type=bigWig',
                                   'name="%s[%s]"' % (args.track_name, strand),
                                   'description="%s[%s]"' % (args.track_desc, strand),
                                   'visibility=full',
                                   'color=%s' % (STRAND_COLOR_MAP[strand]),
                                   'autoScale=on',
                                   'maxHeightPixels=64:64:11',
                                   'bigDataUrl=%s' % (ucsc_url)])
            print >>fileh, track_line
        fileh.close()
    return 0

if __name__ == '__main__':
    sys.exit(main())