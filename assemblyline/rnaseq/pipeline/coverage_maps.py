'''
Created on Dec 10, 2012

@author: mkiyer
'''
import logging
import argparse
import sys
import os
import subprocess

import pysam

import assemblyline.rnaseq.lib.picard as picard
from assemblyline.rnaseq.lib.config import STRAND_SUFFIX_DICT
from assemblyline.rnaseq.lib.base import check_executable
from assemblyline.rnaseq.lib.libtable import FR_UNSTRANDED
from assemblyline.rnaseq.lib.inspect import RnaseqLibraryMetrics

def make_coverage_map(bam_file, bigwig_file, scale_factor, 
                      chrom_sizes_file, tmp_dir):
    # generate bedgraph coverage file
    logging.debug("Generating coverage bedGraph file")
    args = ["bedtools", "genomecov", "-bg", "-split", 
            "-scale", scale_factor,
            "-ibam", bam_file]
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
    logging.debug("Convert bedGraph to bigWig file")        
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

def subdivide_bam_by_strand(bam_file, output_prefix, strand_tag, tmp_dir):
    # setup output files
    strand_bam_files = {}
    strand_bam_fhs = {}
    strand_counts = {'+': 0, '-': 0, '.': 0}
    infh = pysam.Samfile(bam_file, 'rb')
    for strand,suffix in STRAND_SUFFIX_DICT.iteritems():
        filename = os.path.join(tmp_dir, "%s_%s.bam" % 
                                (output_prefix, suffix))
        strand_bam_files[strand] = filename
        strand_bam_fhs[strand] = pysam.Samfile(filename, 'wb', template=infh)
    # iterate through input file
    for r in infh:
        # get strand tag from SAM
        try:
            strand = r.opt(strand_tag)
        except KeyError:
            strand = "."
        # write to appropriate output file
        strand_counts[strand] += 1
        outfh = strand_bam_fhs[strand]
        outfh.write(r)
    infh.close()
    for fh in strand_bam_fhs.itervalues():
        fh.close()
    for strand,count in strand_counts.iteritems():
        logging.debug("Found %d reads on %s strand" % (count, strand))
        if count == 0:
            bam_file = strand_bam_files[strand]
            if os.path.exists(bam_file):
                os.remove(bam_file)
            del strand_bam_files[strand]
    return strand_bam_files

def main():
    # setup logging
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    # command line parsing
    parser = argparse.ArgumentParser()
    parser.add_argument('--scale', dest='scale', action="store_true", default=False)
    parser.add_argument('--strand-tag', dest='strand_tag', default='XS')
    parser.add_argument('--tmp-dir', dest="tmp_dir", default="/tmp")
    parser.add_argument('library_metrics_file')
    parser.add_argument('picard_alignment_metrics_file')
    parser.add_argument('chrom_sizes_file')
    parser.add_argument('bam_file')
    parser.add_argument('bigwig_prefix')
    args = parser.parse_args()
    if not check_executable("bedGraphToBigWig"):
        parser.error("'bedGraphToBigWig' executable not found in PATH")
    if not check_executable("bedtools"):
        parser.error("'bedtools' executable not found in PATH")
    # determine normalization factor from picard metrics
    scale_factor = 1.0 
    if args.scale:
        logging.debug("Reading library size from picard alignment metrics file")
        align_metrics = picard.AlignmentSummaryMetrics(args.picard_alignment_metrics_file)
        library_size = align_metrics.get_total_reads()
        logging.debug("Found library size: %d" % (library_size))
        if library_size > 0:
            scale_factor = 1.0e6 / float(library_size)
    # predict library type
    logging.debug("Getting library type")
    obj = RnaseqLibraryMetrics.from_file(args.library_metrics_file)
    predicted_library_type = obj.predict_library_type()
    logging.debug("Predicted library type: %s" % (predicted_library_type))
    if predicted_library_type == FR_UNSTRANDED:
        # unstranded so make a single coverage map
        logging.debug("Generating unstranded coverage map")
        bigwig_file = args.bigwig_prefix + ".bw"
        retcode = make_coverage_map(args.bam_file, 
                                    bigwig_file, 
                                    scale_factor, 
                                    args.chrom_sizes_file, 
                                    args.tmp_dir)
    else:
        logging.debug("Splitting BAM file by strand")
        prefix = os.path.splitext(os.path.basename(args.bam_file))[0]
        strand_bam_files = subdivide_bam_by_strand(args.bam_file, 
                                                   prefix, 
                                                   args.strand_tag, 
                                                   args.tmp_dir)
        logging.debug("Generating strand-specific coverage maps")
        for strand,bam_file in strand_bam_files.iteritems():
            bigwig_file = "%s_%s.bw" % (args.bigwig_prefix, STRAND_SUFFIX_DICT[strand])
            retcode = make_coverage_map(bam_file, bigwig_file, 
                                        scale_factor, 
                                        args.chrom_sizes_file, 
                                        args.tmp_dir)
            if retcode != 0:
                logging.error("Error")
                if os.path.exists(bigwig_file):
                    os.remove(bigwig_file)
                return 1
        for f in strand_bam_files.itervalues():
            if os.path.exists(f):
                os.remove(f)
    return 0

if __name__ == '__main__':
    sys.exit(main())