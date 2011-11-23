'''
Created on Oct 26, 2011

@author: mkiyer
'''
import sys
import os
import logging
import argparse
import collections
import subprocess
from multiprocessing import Pool

import pysam

from assemblyline.lib.seq import DNA_reverse_complement
from assemblyline.lib.batch_sort import batch_sort
from assemblyline.lib.base import up_to_date
from assemblyline.lib.juncs import parse_raw_juncs

# R script to call for classifying junctions
import assemblyline
_module_dir = assemblyline.__path__[0]
R_SCRIPT = os.path.join(_module_dir, "lib", "classify_juncs.R")

# splice motifs considered "common"
PREVALENT_SPLICE_MOTIFS = set(["GTAG", "GCAG"])
# number of reads considered "enough" for nominating junctions  
ABS_READ_THRESHOLD = 2

class LaneInfo(object):
    __slots__ = ('cohort', 'patient', 'sample', 'lib', 'lane', 'qc',
                 'use_juncs', 'use_transcripts',  'valid', 
                 'aligned_reads', 'read_length', 'tophat_juncs_file')    
    
    @staticmethod
    def from_fields(fields):
        l = LaneInfo()
        l.cohort = fields[0]
        l.patient = fields[1]
        l.sample = fields[2]
        l.lib = fields[3]
        l.lane = fields[4]
        l.qc = fields[5]
        l.use_juncs = int(fields[6])
        l.use_transcripts = int(fields[7])        
        l.valid = int(fields[8])
        if fields[9] == "NA":
            l.aligned_reads = None
        else:
            l.aligned_reads = int(float(fields[9]))
        l.read_length = int(fields[10])
        if fields[11] == "NA":
            l.tophat_juncs_file = None
        else:
            l.tophat_juncs_file = fields[11]
        return l

    def is_valid(self):
        if self.qc == "FAIL":
            return False
        if self.aligned_reads is None:
            return False
        if self.use_juncs and (self.tophat_juncs_file is not None):
            return True
        return False

    @staticmethod
    def from_file(filename):
        fh = open(filename)
        # skip header
        fh.next()
        for line in fh:
            fields = line.strip().split('\t')
            yield LaneInfo.from_fields(fields)
        fh.close()

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

class AnnotatedJunctionData(object):
    __slots__ = ('chrom', 'start', 'end', 'strand', 'enough_reads',
                 'best_cov', 'best_overhang', 'motif', 'recur', 'category', 'annotated')
    @staticmethod
    def from_fields(fields):
        j = AnnotatedJunctionData()
        j.chrom = fields[0]
        j.start = int(fields[1])
        j.end = int(fields[2])
        j.strand = fields[3]
        j.enough_reads = int(fields[4])
        j.best_cov = float(fields[5])
        j.best_overhang = int(fields[6])
        j.motif = int(fields[7])
        j.recur = int(fields[8])
        j.category = fields[9]
        j.annotated = int(fields[10])
        return j

    def to_line(self):
        return '\t'.join(map(str, [self.chrom, self.start, self.end, self.strand,
                                   self.enough_reads, self.best_cov, self.best_overhang,
                                   self.motif, self.recur, self.category, self.annotated]))

class SummaryJunctionData(object):
    __slots__ = ('chrom', 'start', 'end', 'strand', 'motif', 'recur', 
                 'category', 'annotated', 'lanes_true', 'lanes_false')

    @staticmethod
    def from_fields(fields):
        j = SummaryJunctionData()
        j.chrom = fields[0]
        j.start = int(fields[1])
        j.end = int(fields[2])
        j.strand = fields[3]
        j.motif = int(fields[4])
        j.recur = int(fields[5])
        j.category = fields[6]
        j.annotated = int(fields[7])
        j.lanes_true = int(fields[8])
        j.lanes_false = int(fields[9])

    @staticmethod
    def from_file(filename):
        fh = open(filename)
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.strip().split('\t')
            yield LaneInfo.from_fields(fields)
        fh.close()

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

def aggregate_juncs_data(laneinfos, junc_data_file,
                         tmp_dir="/tmp"):
    """
    aggregates junctions from all samples from the 'sample_table_file', 
    computes coverage, and outputs results sorted by genomic position
    """
    outfh = open(junc_data_file, "w")
    for laneinfo in laneinfos:
        logging.debug("\tcollecting junctions from %s|%s|%s|%s|%s" % 
                      (laneinfo.cohort, laneinfo.patient, laneinfo.sample, 
                       laneinfo.lib, laneinfo.lane))
        j = JunctionData()
        j.cohort = laneinfo.cohort
        j.patient = laneinfo.patient
        j.sample = laneinfo.sample
        j.lib = laneinfo.lib
        j.lane = laneinfo.lane
        for chrom, left, right, strand, score, overhang in \
            parse_tophat_junctions_bed(open(laneinfo.tophat_juncs_file)):
            j.chrom = chrom
            j.start = left
            j.end = right
            j.strand = strand
            j.overhang = overhang
            j.reads = score
            # compute coverage
            j.cov = float(score)*(1.0e3/laneinfo.read_length)*(1.0e6/laneinfo.aligned_reads)
            # output information to master file
            print >>outfh, j.to_line()
    outfh.close()


def get_junc_category(chrom, start, end, strand, known_juncs, 
                      known_left_splice_sites, known_right_splice_sites):
    is_known = int((chrom, start, end, strand) in known_juncs)
    is_left_known = (chrom, start, strand) in known_left_splice_sites
    is_right_known = (chrom, end, strand) in known_right_splice_sites
    if is_known:
        return "aa"
    elif (is_left_known and is_right_known):
        return "kk"
    elif (is_left_known or is_right_known):
        return "kn"
    else:
        return "nn"

def make_annot_junc(chrom, start, end, strand, patients, 
                    best_reads, best_cov, best_overhang, 
                    ref_fa, known_juncs, known_left_splice_sites, 
                    known_right_splice_sites):
    # check whether genomic splice motif is common
    motif = get_splice_motif(chrom, start, end, strand, ref_fa)
    is_prevalent_motif = int(motif in PREVALENT_SPLICE_MOTIFS)
    # get junction category
    is_known = int((chrom, start, end, strand) in known_juncs)
    category = get_junc_category(chrom, start, end, strand,
                                 known_juncs, known_left_splice_sites,
                                 known_right_splice_sites)
    # get recurrence of junction
    patient_recur = len(patients)
    # see whether junction occurs in enough reads
    enough_reads = int(best_reads >= ABS_READ_THRESHOLD)
    # print
    j = AnnotatedJunctionData()
    j.chrom = chrom
    j.start = start
    j.end = end
    j.strand = strand
    j.enough_reads = enough_reads
    j.best_cov = best_cov
    j.best_overhang = best_overhang
    j.motif = is_prevalent_motif
    j.recur = patient_recur
    j.category = category
    j.annotated = is_known
    return j

def condense_annotate_juncs(junc_data_file, output_file, ref_fasta_file, 
                            known_raw_juncs_file):
    """
    given a sorted aggregated junctions file supplied by the collect_juncs_data
    function, annotates each junction with number of reads, coverage, overhang, and
    recurrence across multiple samples
    
    creates two output files: one for "known junctions" (training set) and one for
    unannotated junctions (test set)
    """ 
    # read known junctions
    logging.debug("\treading list of known junctions")
    known_juncs = set()
    known_left_splice_sites = set()
    known_right_splice_sites = set()    
    for chrom, start, end, strand in parse_raw_juncs(open(known_raw_juncs_file)):
        known_juncs.add((chrom, start, end, strand))
        known_left_splice_sites.add((chrom, start, strand))
        known_right_splice_sites.add((chrom, end, strand))        
    # prepare output files    
    outfh = open(output_file, "w")
    header_fields = ["#chrom", "start", "end", "strand", "enough_reads", 
                     "best_cov", "best_overhang", "motif", "recur", "category", "annotated"]   
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
                aj = make_annot_junc(prev[0], prev[1], prev[2], prev[3],
                                     patients, best_reads, best_cov, best_overhang, 
                                     ref_fa, known_juncs, known_left_splice_sites, 
                                     known_right_splice_sites)                    
                print >>outfh, aj.to_line()
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
        aj = make_annot_junc(j.chrom, j.start, j.end, j.strand,
                             patients, best_reads, best_cov, best_overhang, 
                             ref_fa, known_juncs, known_left_splice_sites, 
                             known_right_splice_sites)
        print >>outfh, aj.to_line()
    ref_fa.close()
    outfh.close()

def split_juncs_by_category(consensus_junc_file, outfiles):
    # open output files
    outfhs = {}
    counts = collections.defaultdict(lambda: 0)
    fh = open(consensus_junc_file)
    header = fh.next().strip()    
    for category in ("kk", "kn", "nn"):
        outfhs[category] = open(outfiles[category], "w")
        # write header
        print >>outfhs[category], header
    for line in fh:
        fields = line.strip().split('\t')
        aj = AnnotatedJunctionData.from_fields(fields)
        if aj.category == "aa":
            for outfh in outfhs.itervalues():
                print >>outfh, aj.to_line()
        else:
            counts[aj.category] += 1
            print >>outfhs[aj.category], aj.to_line()
    for outfh in outfhs.itervalues():
        outfh.close()
    fh.close()
    return outfiles, counts

def get_junc_ctree_files(laneinfo, output_dir):
    prefix = "%s_%s" % (laneinfo.sample, laneinfo.lane)
    outprefix = os.path.join(output_dir, prefix)
    juncfile = outprefix + ".juncs"
    outfile = outprefix + ".classify.juncs"
    return juncfile, outprefix, outfile

def run_junc_classifier(args):
    infile, outprefix = args
    logfile = outprefix + ".classify.log"
    logging.debug("\trunning classifier on file %s" % (infile))
    logfh = open(logfile, "w")
    retcode = subprocess.call(["Rscript", R_SCRIPT, infile, outprefix], stdout=logfh, stderr=logfh)
    logfh.close()
    if retcode != 0:
        logging.error("\tclassification failed with error code %d" % (retcode))
    return (retcode, outprefix)

def classify_juncs(laneinfos, consensus_juncs_file, output_dir,
                   num_processors=1):
    # load master juncs into memory
    logging.debug("Reading consensus junctions into memory")
    annot_juncs = {}
    # skip header
    fh = open(consensus_juncs_file)
    fh.next().strip().split('\t')
    for line in fh:
        fields = line.strip().split('\t')
        aj = AnnotatedJunctionData.from_fields(fields)
        annot_juncs[(aj.chrom, aj.start, aj.end, aj.strand)] = aj
    fh.close()
    # write lane classification files one by one
    tasks = []
    for laneinfo in laneinfos:
        # input file
        juncfile, outprefix, outfile = get_junc_ctree_files(laneinfo, output_dir)
        msg = ("\twriting classification file for cohort=%s pt=%s sample=%s lib=%s lane=%s" % 
               (laneinfo.cohort, laneinfo.patient, laneinfo.sample, 
                laneinfo.lib, laneinfo.lane))
        if up_to_date(juncfile, laneinfo.tophat_juncs_file):
            logging.debug("[SKIPPED] %s" % (msg))
        else:
            logging.debug(msg)
            juncfh = open(juncfile, "w")
            print >>juncfh, '\t'.join(["chrom", "start", "end", "strand",
                                       "myreads", "myoverhang", "enough_reads", 
                                       "best_cov", "best_overhang",
                                       "motif", "recur", "category", "annotated"])
            for chrom, left, right, strand, score, overhang in \
                parse_tophat_junctions_bed(open(laneinfo.tophat_juncs_file)):
                # lookup annotated junc
                if (chrom, left, right, strand) not in annot_juncs:
                    continue
                aj = annot_juncs[(chrom, left, right, strand)]
                # output junc            
                print >>juncfh, '\t'.join(map(str, [chrom, left, right, strand, 
                                                    score, overhang, aj.enough_reads,
                                                    aj.best_cov, aj.best_overhang,
                                                    aj.motif, aj.recur, aj.category, 
                                                    aj.annotated]))
            juncfh.close()
        msg = ("\tadding classification task for cohort=%s pt=%s sample=%s lib=%s lane=%s" % 
               (laneinfo.cohort, laneinfo.patient, laneinfo.sample, 
                laneinfo.lib, laneinfo.lane))
        if up_to_date(outfile, juncfile):
            logging.debug("[SKIPPED] %s" % (msg))
        else:
            logging.debug(msg)
            # add classification task
            tasks.append((juncfile, outprefix))
    # use multiprocessing to parallelize classification
    pool = Pool(processes=num_processors)
    result_iter = pool.imap_unordered(run_junc_classifier, tasks)
    for retcode,outprefix in result_iter:
        if retcode == 0:
            logging.debug("\tsuccessfully classified %s" % (outprefix))
        else:
            logging.debug("\terror when classifying %s" % (outprefix))


def summarize_classification_results(laneinfos, predicted_junc_file, junc_ctree_dir):
    # summarize junction classification results
    juncs = {}
    for laneinfo in laneinfos:
        # construct filenames of output
        juncfile, outprefix, outfile = get_junc_ctree_files(laneinfo, junc_ctree_dir)
        cutfile = outprefix + ".classify.cutoff"
        # init stats
        tp, fp, tn, fn = 0, 0, 0, 0
        # read cutoff
        cutoff = float(open(cutfile).next().strip())    
        # read classification tree output
        fh = open(outfile)
        header_fields = fh.next().strip().split('\t')
        annotated_index = header_fields.index("annotated")
        motif_index = header_fields.index("motif")
        recur_index = header_fields.index("recur")
        result_index = header_fields.index("result")
        for line in fh:
            # read line
            fields = line.strip().split('\t')
            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            strand = fields[3]
            motif = int(fields[motif_index])
            recur = int(fields[recur_index])
            annotated = int(fields[annotated_index])
            prob = float(fields[result_index])
            # insert into dictionary
            key = (chrom, start, end, strand)
            if (chrom, start, end, strand) in juncs:
                j = juncs[key]
            else:
                j = SummaryJunctionData()
                j.chrom = chrom
                j.start = start
                j.end = end
                j.strand = strand
                j.motif = motif
                j.recur = recur
                j.annotated = annotated
                j.lanes_true = 0
                j.lanes_false = 0
                juncs[key] = j
            # keep track of statistics
            if prob >= cutoff:
                j.lanes_true += 1
                if annotated:
                    tp += 1
                else:
                    fp += 1
            else:
                j.lanes_false += 1
                if annotated:
                    fn += 1
                else:
                    tn += 1
        # output lane results
        prec = tp / float(tp + fp)
        rec = tp / float(tp + fn)
        spec = (tn / float(tn + fp))
        sens = (tp / float(tp + fn))
        acc = (tp + tn) / float(tp + fp + tn + fn)
        if tp == 0:
            F = 0
        else:
            F = 2 * (prec * rec) / (prec + rec)
        logging.info("Lane=%s Cutoff=%f TP=%d FN=%d FP=%d TN=%d Prec=%f Recall/Sens=%f Spec=%f F=%f Acc=%f" % 
                     (laneinfo.lane, cutoff, tp, fn, fp, tn, prec, sens, spec, F, acc))
    # output junction summary data
    tp, fp, tn, fn = 0, 0, 0, 0
    fh = open(predicted_junc_file, "w")
    print >>fh, '\t'.join(["#chrom", "start", "end", "strand", "motif", 
                           "recur", "annotated", "lanes_true", 
                           "lanes_false"])
    for k in sorted(juncs):
        j = juncs[k]
        # keep track of statistics
        if j.lanes_true > 0:
            if j.annotated:
                tp += 1
            else:
                fp += 1
        elif j.annotated:
            fn += 1
        else:
            tn += 1
        # output final summary report
        print >>fh, '\t'.join(map(str, [j.chrom, j.start, j.end, j.strand, 
                                        j.motif, j.recur, j.annotated, 
                                        j.lanes_true, j.lanes_false]))
    fh.close()

    prec = tp / float(tp + fp)
    rec = tp / float(tp + fn)
    sens = rec
    spec = (tn / float(tn + fp))
    acc = (tp + tn) / float(tp + fp + tn + fn)
    if tp == 0:
        F = 0
    else:
        F = 2 * (prec * rec) / (prec + rec)    
    logging.info("TP=%d FN=%d FP=%d TN=%d" % (tp, fn, fp, tn))
    logging.info("Precision=%f Recall/Sensitivity=%f" % (prec, rec))
    logging.info("Specificity=%f" % (spec))
    logging.info("F=%f" % F)
    logging.info("Accuracy=%f" % (acc))

def main():
    #
    # Command line parsing
    # 
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--output-dir', dest="output_dir", default=".")
    parser.add_argument('-p', '--num-processors', dest="num_processors", type=int, default=1)
    parser.add_argument('reference_genome')
    parser.add_argument('known_raw_juncs_file')
    parser.add_argument('sample_table')
    args = parser.parse_args()
    #
    # Setup output locations
    #
    if not os.path.exists(args.output_dir):
        logging.info("Creating output directory: %s" % (args.output_dir))
        os.makedirs(args.output_dir)
    tmp_dir = os.path.join(args.output_dir, "tmp")
    if not os.path.exists(tmp_dir):
        logging.info("Creating tmp directory: %s" % (tmp_dir))
        os.makedirs(tmp_dir)
    # setup intermediate file names
    junc_data_file = os.path.join(tmp_dir, "aggregated_juncs.txt")
    sorted_junc_data_file = os.path.join(tmp_dir, "aggregated_juncs.srt.txt")
    consensus_junc_file = os.path.join(args.output_dir, "annotated_consensus_juncs.txt")
    #
    # Parse sample table
    #
    logging.info("Parsing sample table")
    laneinfos = []
    for laneinfo in LaneInfo.from_file(args.sample_table):
        if not laneinfo.is_valid():
            logging.warning("\tskipping %s|%s|%s|%s|%s" % 
                            (laneinfo.cohort, laneinfo.patient, laneinfo.sample, 
                             laneinfo.lib, laneinfo.lane))
            continue
        laneinfos.append(laneinfo)
    #
    # Aggregate junction information from all runs
    #
    msg = "Aggregating junction data"
    if up_to_date(junc_data_file, args.sample_table):
        logging.info("[SKIPPED] %s" % (msg))
    else:    
        logging.info(msg)
        aggregate_juncs_data(laneinfos,
                             junc_data_file,
                             tmp_dir=tmp_dir)
    #
    # Sort aggregated junction list
    #
    msg = "Sorting junctions by genomic position"
    if up_to_date(sorted_junc_data_file, junc_data_file):
        logging.info("[SKIPPED] %s" % (msg))
    else:
        logging.info(msg)    
        sort_juncs_by_pos(junc_data_file, sorted_junc_data_file, tmp_dir)
    #
    # Annotate and condense junctions
    #
    msg = "Condensing and annotating junctions"
    if up_to_date(consensus_junc_file, sorted_junc_data_file):
        logging.info("[SKIPPED] %s" % msg)
    else:
        logging.info(msg)
        condense_annotate_juncs(sorted_junc_data_file, consensus_junc_file, 
                                args.reference_genome, args.known_raw_juncs_file)
    #
    # Split juncs by category
    #
    msg = "Splitting junctions into categories"
    consensus_category_files = {}
    for category in ("kk", "kn", "nn"):
        consensus_category_files[category] = os.path.join(args.output_dir, "annotated_consensus_juncs.%s.txt" % (category))
    skip = True    
    for filename in consensus_category_files.itervalues():    
        if not up_to_date(filename, consensus_junc_file):
            skip = False
            break
    if skip:
        logging.info("[SKIPPED] %s" % (msg))
    else:
        logging.info(msg)
        split_juncs_by_category(consensus_junc_file, consensus_category_files)
    #
    # Classify junctions per category
    # 
    for category,consensus_category_file in consensus_category_files.iteritems():
        ctree_dir = os.path.join(args.output_dir, "ctrees_%s" % (category))
        if not os.path.exists(ctree_dir):
            logging.info("Creating result directory: %s" % (ctree_dir))
            os.makedirs(ctree_dir)
        # Classify junctions
        classify_juncs(laneinfos, consensus_category_file, ctree_dir, args.num_processors)
        # Summarize classification results
        predicted_junc_file = os.path.join(args.output_dir, "predicted_juncs.%s.txt" % (category))
        msg = "Summarizing classification results"
        skip = True
        for laneinfo in laneinfos:
            juncfile, outprefix, outfile = get_junc_ctree_files(laneinfo, ctree_dir)
            if not up_to_date(predicted_junc_file, outfile):
                skip = False
                break
        if skip:
            logging.info("[SKIPPED] %s" % (msg))
        else:
            logging.info(msg)
            summarize_classification_results(laneinfos, predicted_junc_file, ctree_dir)

if __name__ == '__main__':
    sys.exit(main())
