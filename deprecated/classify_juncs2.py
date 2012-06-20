'''
Created on Nov 11, 2011

@author: mkiyer
'''
import sys
import os
import logging
import argparse
import pysam
import subprocess

from assemblyline.lib.seq import DNA_reverse_complement
from assemblyline.lib.batch_sort import batch_sort

# R script to call for classifying junctions
import assemblyline
_module_dir = assemblyline.__path__[0]
R_SCRIPT = os.path.join(_module_dir, "lib", "classify_juncs.R")

# splice motifs considered "common"
PREVALENT_SPLICE_MOTIFS = set(["GTAG", "GCAG"])
# number of reads considered "enough" for nominating junctions  
ABS_READ_THRESHOLD = 2

class JunctionData(object):
    """
    holds junction data from a single lane
    """
    __slots__ = ('chrom', 'start', 'end', 'strand', 'reads', 'cov', 
                 'overhang', 'cohort', 'patient', 'sample', 'lib', 'lane')
    @staticmethod
    def from_fields(fields):
        j = JunctionData()
        j.chrom = fields[0]
        j.start = int(fields[1])
        j.end = int(fields[2])
        j.strand = fields[3]
        j.reads = int(fields[4])
        j.cov = float(fields[5])
        j.overhang = int(fields[6])
        j.cohort = fields[7]
        j.patient = fields[8]
        j.sample = fields[9]
        j.lib = fields[10]
        j.lane = fields[11]
        return j
    
    def to_line(self):
        return '\t'.join(map(str, [self.chrom, self.start, self.end, self.strand,
                                   self.reads, self.cov, self.overhang,
                                   self.cohort, self.patient, self.sample, self.lib, self.lane]))


def get_splice_motif(chrom, start, end, strand, ref_fa):
    """
    returns the genomic DNA splicing motif for a given junction
    """
    left = ref_fa.fetch(chrom, start+1, start+3)
    right = ref_fa.fetch(chrom, end-2, end)
    motif = left + right
    if strand == "-":
        motif = DNA_reverse_complement(motif)
    return motif
    
def parse_raw_juncs(line_iter):
    """
    generator to parse lines in the tophat "raw juncs" format
    """
    for line in line_iter:
        fields = line.strip().split('\t')
        chrom = fields[0]
        start = int(fields[1])
        end = int(fields[2])
        strand = fields[3]
        yield (chrom, start, end, strand)
        
def parse_tophat_junctions_bed(line_iter):
    """
    generator to parse tophat junctions bed file and return junction data
    as tuples
    """
    for linenum,line in enumerate(line_iter):
        if not line:
            continue
        if line.startswith("#") or line.startswith("track"):
            continue
        fields = line.strip().split('\t')
        if len(fields) < 12:
            logging.warning("malformed line %d, missing columns" % linenum)
            continue
        chrom = fields[0]
        start = int(fields[1])
        score = int(fields[4])
        strand = fields[5]
        block_sizes = map(int, [x for x in fields[10].split(',') if x])
        block_starts = map(int, [x for x in fields[11].split(',') if x])
        # get smallest overhang
        overhang = min(block_sizes)
        for i in xrange(1,len(block_sizes)):
            left_pos = start + block_starts[i-1] + block_sizes[i-1] - 1
            right_pos = start + block_starts[i]            
            yield chrom, left_pos, right_pos, strand, score, overhang

def sort_juncs_by_pos(input_file, output_file, tmp_dir):
    def sort_junc_func(line):
        fields = line.strip().split('\t', 4)[:4]
        return fields[0], int(fields[1]), int(fields[2]), fields[3]
    batch_sort(input=input_file,
               output=output_file,
               key=sort_junc_func,
               buffer_size=(1 << 21),
               tempdirs=[tmp_dir])

def collect_juncs_data(sample_table_file,
                       junc_data_file,
                       tmp_dir="/tmp"):
    """
    aggregates junctions from all samples from the 'sample_table_file', 
    computes coverage, and outputs results sorted by genomic position
    """
    # read sample table
    sample_table_fh = open(sample_table_file)
    header_fields = sample_table_fh.next().strip().split('\t')
    junc_bed_index = header_fields.index("tophat_juncs_file")
    aligned_reads_index = header_fields.index("aligned_reads")
    read_length_index = header_fields.index("read_length")    
    outfh = open(junc_data_file, "w")
    for line in sample_table_fh:
        # cohort, patient, sample, library, lane, aligned_reads, read_length, tophat_juncs_file    
        fields = line.strip().split('\t')
        j = JunctionData()
        j.cohort = fields[0]
        j.patient = fields[1]
        j.sample = fields[2]
        j.lib = fields[3]
        j.lane = fields[4]
        junctions_bed_file = fields[junc_bed_index]
        aligned_reads = int(float(fields[aligned_reads_index]))
        read_length = int(fields[read_length_index])
        logging.info("Collecting junctions from %s|%s|%s|%s|%s" % (j.cohort, j.patient, j.sample, j.lib, j.lane))
        for chrom, left, right, strand, score, overhang in parse_tophat_junctions_bed(open(junctions_bed_file)):
            j.chrom = chrom
            j.start = left
            j.end = right
            j.strand = strand
            j.overhang = overhang
            j.reads = score
            # compute coverage
            j.cov = float(score)*(1.0e3/read_length)*(1.0e6/aligned_reads)
            # output information to master file
            print >>outfh, j.to_line()
    sample_table_fh.close()
    outfh.close()

def output_junc_info(chrom, start, end, strand, patients, reads, cov, 
                     overhang, ref_fa, known_juncs):
    # check whether genomic splice motif is common
    motif = get_splice_motif(chrom, start, end, strand, ref_fa)
    is_prevalent_motif = int(motif in PREVALENT_SPLICE_MOTIFS)
    # check whether junc is known
    is_known = int((chrom, start, end, strand) in known_juncs)
    # get recurrence of junction
    patient_recur = len(patients)
    # see whether junction occurs in enough reads
    enough_reads = int(reads >= ABS_READ_THRESHOLD)
    # print
    fields = '\t'.join(map(str, [chrom, start, end, strand,
                                 enough_reads, cov, overhang,
                                 is_prevalent_motif, patient_recur, 
                                 is_known]))
    return fields

def annotate_juncs(junc_data_file, output_file, ref_fasta_file, 
                   known_raw_juncs_file):
    """
    given a sorted aggregated junctions file supplied by the collect_juncs_data
    function, annotates each junction with number of reads, coverage, overhang, and
    recurrence across multiple samples
    """ 
    # read known junctions
    logging.debug("\treading list of known junctions")
    known_juncs = set()
    for chrom, start, end, strand in parse_raw_juncs(open(known_raw_juncs_file)):
        known_juncs.add((chrom, start, end, strand))
    # prepare output files    
    outfh = open(output_file, "w")
    header_fields = ["chrom", "start", "end", "strand", "best_reads", 
                     "best_cov", "best_overhang", "motif", "recur", "annotated"]   
    print >>outfh, '\t'.join(header_fields)
    # read junction data
    patients = set()
    best_reads = 0
    best_cov = 0
    best_overhang = 0
    prev = None
    ref_fa = pysam.Fastafile(ref_fasta_file)
    logging.debug("\tcondensing and annotating junctions")
    for line in open(junc_data_file):
        j = JunctionData.from_fields(line.strip().split('\t'))
        if (prev != (j.chrom, j.start, j.end, j.strand)):
            if prev is not None:
                print >>outfh, output_junc_info(prev[0], prev[1], prev[2], prev[3],
                                                patients, best_reads, best_cov, best_overhang, 
                                                ref_fa, known_juncs)
            # reset statistics
            prev = (j.chrom, j.start, j.end, j.strand)
            patients = set()
            best_reads = 0
            best_cov = 0
            best_overhang = 0
        # update stats
        patients.add(j.patient)
        best_reads = max(best_reads, j.reads)
        best_cov = max(best_cov, j.cov)
        best_overhang = max(best_overhang, j.overhang) 
    if len(patients) > 0:
        output_junc_info(j.chrom, j.start, j.end, j.strand,
                         patients, best_reads, best_cov, best_overhang, 
                         ref_fa, known_juncs)
    ref_fa.close()
    outfh.close()

def classify_juncs(infile, output_dir):
    prefix = os.path.splitext(os.path.basename(infile))[0]
    outfile = os.path.join(output_dir, prefix + ".classify.txt")
    cutfile = os.path.join(output_dir, prefix + ".classify.cutoff")
    epsfile = os.path.join(output_dir, prefix + ".classify.eps")
    logfile = os.path.join(output_dir, prefix + ".classify.log")
    logfh = open(logfile, "w")
    retcode = subprocess.call(["Rscript", R_SCRIPT, infile, outfile, 
                               cutfile, epsfile], stdout=logfh, 
                               stderr=logfh)
    logfh.close()
    if retcode != 0:
        logging.error("\tclassification failed with error code %d" % (retcode))
    # read cutoff
    cutoff = float(open(cutfile).next().strip())    
    logging.debug("\tusing cutoff probability %f" % (cutoff))
    # write output
    tp, fp, tn, fn = 0, 0, 0, 0
    fh = open(outfile)
    header_fields = fh.next().strip().split('\t')
    annotated_index = header_fields.index("annotated")
    result_index = header_fields.index("result")
    goodfh = open(os.path.join(output_dir, "good_juncs.txt"), "w")
    badfh = open(os.path.join(output_dir, "bad_juncs.txt"), "w")
    for line in fh:
        fields = line.strip().split('\t')
        prob = float(fields[result_index])
        annotated = int(fields[annotated_index])
        # keep track of statistics
        if annotated and (prob >= cutoff):
            tp += 1
        elif annotated:
            fn += 1
        elif (prob >= cutoff):
            fp += 1
        else:
            tn += 1
        if annotated or (prob >= cutoff):
            print >>goodfh, '\t'.join(fields)
        else:
            print >>badfh, '\t'.join(fields)
    spec = (tn / float(tn + fp))
    sens = (tp / float(tp + fn))
    acc = 2 * (spec * sens) / (sens + spec)
    logging.info("TP=%d FN=%d FP=%d TN=%d" % (tp, fn, fp, tn))
    logging.info("Sensitivity: %f" % (spec))
    logging.info("Specificity: %f" % (sens))
    logging.info("Accuracy: %f" % (acc))

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--output-dir', dest="output_dir", default=".")
    parser.add_argument('--tmp-dir', dest="tmp_dir", default="/tmp")
    parser.add_argument('reference_genome')
    parser.add_argument('known_raw_juncs_file')
    parser.add_argument('sample_table')
    args = parser.parse_args()
    junc_data_file = os.path.join(args.output_dir, "aggregated_juncs.txt")
    sorted_junc_data_file = os.path.join(args.output_dir, "aggregated_juncs.srt.txt")
    consensus_junc_file = os.path.join(args.output_dir, "consensus_juncs.srt.txt")
    logging.info("Aggregating junction data...")
    collect_juncs_data(args.sample_table,
                       junc_data_file,
                       tmp_dir=args.tmp_dir)
    logging.info("Sorting junctions by genomic position")
    sort_juncs_by_pos(junc_data_file, sorted_junc_data_file, args.tmp_dir)
    logging.info("Annotating junctions...")
    annotate_juncs(sorted_junc_data_file, consensus_junc_file, args.sample_table, 
                   args.reference_genome, args.known_raw_juncs_file)
    logging.info("Classifying unannotated junctions...")
    classify_juncs(consensus_junc_file, args.output_dir)


if __name__ == '__main__':
    sys.exit(main())
