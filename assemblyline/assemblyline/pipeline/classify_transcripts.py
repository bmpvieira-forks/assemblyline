'''
Created on Nov 30, 2011

@author: mkiyer
'''
import logging
import argparse
import os
import collections
import subprocess
import shutil
import heapq
import sys
import numpy as np

from assemblyline.lib.transcript import NO_STRAND, POS_STRAND, NEG_STRAND
from assemblyline.lib.transcript_parser import parse_gtf
from assemblyline.lib.sampletable import LibInfo
from assemblyline.lib.bx.cluster import ClusterTree
from assemblyline.lib.bx.intersection import Interval, IntervalTree
from assemblyline.lib.batch_sort import batch_sort
from assemblyline.lib import bed
from assemblyline.lib import gtf

# constant transcript attribute names
DENSITY = 'd'
RECUR = 'r'
SENSE_EXON_FRAC = 'sf'
ANTISENSE_EXON_FRAC = 'asf'
INTRON_FRAC = 'if'

# overlap necessary to consider a transcript "annotated"
MIN_SENSE_EXON_FRAC = 0.9

# R script to call for classifying transcripts
import assemblyline
_module_dir = assemblyline.__path__[0]
R_SCRIPT = os.path.join(_module_dir, "lib", "classify_transcripts.R")

def build_known_transcript_trees(bed_iter):
    exon_trees = {"+": collections.defaultdict(lambda: IntervalTree()),
                  "-": collections.defaultdict(lambda: IntervalTree())}
    intron_trees = collections.defaultdict(lambda: IntervalTree())
    # cluster gene boundaries
    locus_cluster_trees = collections.defaultdict(lambda: ClusterTree(0, 1))
    genes = []
    for g in bed.BEDFeature.parse(bed_iter):
        locus_cluster_trees[g.chrom].insert(g.tx_start, g.tx_end, len(genes))
        genes.append(g)
    for chrom, locus_cluster_tree in locus_cluster_trees.iteritems():
        loci = []
        for locus_start, locus_end, indexes in locus_cluster_tree.getregions():        
            loci.append((locus_start, locus_end))
            # cluster exons within gene "locus"
            exon_cluster_trees = {"+": ClusterTree(0,1),
                                  "-": ClusterTree(0,1)}
            intron_cluster_tree = ClusterTree(0,1)
            chrom = genes[indexes[0]].chrom
            for index in indexes:
                g = genes[index]
                for start, end in g.exons:
                    exon_cluster_trees[g.strand].insert(start, end, 0)
                    intron_cluster_tree.insert(start, end, 0)
            # now get clusters and insert into exon trees
            for strand, exon_cluster_tree in exon_cluster_trees.iteritems():            
                for exon_start, exon_end, indexes in exon_cluster_tree.getregions():
                    exon_trees[strand][chrom].insert_interval(Interval(exon_start, exon_end))
            # insert purely intronic regions into intron trees
            exon_boundaries = [(s,e) for s,e,i in intron_cluster_tree.getregions()]
            for e1,e2 in zip(exon_boundaries[:-1], exon_boundaries[1:]):
                intron_trees[chrom].insert_interval(Interval(e1[1], e2[0]))
    return exon_trees, intron_trees

def get_overlap_fraction(tx, trees):
    # Determine percent overlap with annotated genes
    overlap_fraction = 0.0
    if tx.chrom in trees:
        tree = trees[tx.chrom]            
        overlap_bases = 0
        tx_length = 0    
        for e in tx.exons:  
            astart = e.start
            aend = e.end              
            bintervals = tree.find(astart, aend)
            max_exon_overlap_bases = 0
            for interval in bintervals:            
                bstart, bend = interval.start, interval.end
                overlap_start = astart if bstart < astart else bstart
                overlap_end = bend if bend < aend else aend
                max_exon_overlap_bases = max(max_exon_overlap_bases, (overlap_end - overlap_start))
            overlap_bases += max_exon_overlap_bases
            tx_length += (aend - astart)            
        overlap_fraction = overlap_bases / float(tx_length)
    return overlap_fraction

def annotate_overlap(transcripts, exon_trees, intron_trees):
    for tx in transcripts:
        pos_exonic_frac = get_overlap_fraction(tx, exon_trees["+"])
        neg_exonic_frac = get_overlap_fraction(tx, exon_trees["-"])
        intronic_frac = get_overlap_fraction(tx, intron_trees)
        if tx.strand == POS_STRAND:
            sense_frac, antisense_frac = pos_exonic_frac, neg_exonic_frac
        elif tx.strand == NEG_STRAND:
            sense_frac, antisense_frac = neg_exonic_frac, pos_exonic_frac
        else:
            sense_frac = max(pos_exonic_frac, neg_exonic_frac)
            antisense_frac = 0.0
        tx.attrs[SENSE_EXON_FRAC] = sense_frac
        tx.attrs[ANTISENSE_EXON_FRAC] = antisense_frac
        tx.attrs[INTRON_FRAC] = intronic_frac        

def annotate_recurrence_and_density(transcripts, libs_per_sample):
    # group transcripts by sample
    sample_tx_dict = collections.defaultdict(lambda: [])
    for tx in transcripts:
        sample_tx_dict[tx.sample].append(tx)
    # create locus arrays
    locus_start = transcripts[0].start
    locus_end = max(tx.end for tx in transcripts)
    recur_arr = np.zeros((3,(locus_end - locus_start)), dtype=int)
    mass_arr = np.zeros((3,(locus_end - locus_start)), dtype=float)
    # process transcripts from each sample
    for sample, sample_transcripts in sample_tx_dict.iteritems():        
        sample_recur_arr = np.zeros((3,(locus_end - locus_start)), dtype=float)
        sample_mass_arr = np.zeros((3,(locus_end - locus_start)), dtype=float)
        for tx in sample_transcripts:        
            for e in tx.exons:
                sample_recur_arr[tx.strand, (e.start - locus_start):(e.end - locus_start)] = 1
                sample_mass_arr[tx.strand, (e.start - locus_start):(e.end - locus_start)] += tx.mass
        # average the mass arrays to account for multiple libraries per sample
        sample_mass_arr /= libs_per_sample[sample]
        # add sample arrays to compendia arrays
        recur_arr += sample_recur_arr
        mass_arr += sample_mass_arr
    # compute the fraction of positive strand density at each base
    total_stranded_mass_arr = mass_arr[POS_STRAND] + mass_arr[NEG_STRAND]
    nzinds = np.flatnonzero(total_stranded_mass_arr)
    pos_strand_frac_arr = np.ones((locus_end - locus_start), dtype=float)
    pos_strand_frac_arr[nzinds] = mass_arr[POS_STRAND, nzinds] / total_stranded_mass_arr[nzinds]                
    # reallocate unstranded mass proportionately
    strand_frac_arrs = (pos_strand_frac_arr, 1.0 - pos_strand_frac_arr)
    for strand in (POS_STRAND, NEG_STRAND):
        mass_arr[strand] += strand_frac_arrs[strand] * mass_arr[NO_STRAND]
        recur_arr[strand] = np.max(np.vstack((recur_arr[strand], strand_frac_arrs[strand] * recur_arr[NO_STRAND])), axis=0)
    # compute recurrence/density of individual transcripts using the arrays
    for tx in transcripts:
        # sum mass on all strands over the exonic intervals of the transcript
        mass_sums = np.zeros(3, dtype=float)
        recur_sums = np.zeros(3, dtype=float)
        for e in tx.exons:           
            mass_sums += np.sum(mass_arr[:,(e.start - locus_start):(e.end - locus_start)], axis=1)
            recur_sums += np.sum(recur_arr[:,(e.start - locus_start):(e.end - locus_start)], axis=1)
        # assign most likely strand to unstranded transcripts
        if tx.strand == NO_STRAND:
            assigned_strand = np.argmax(mass_sums[0:2])
        else:
            assigned_strand = tx.strand
        # TODO: remove
        assert assigned_strand != NO_STRAND
        # calculate density (frac*mass/length)
        density = tx.frac * mass_sums[assigned_strand] / float(tx.length)
        recur = recur_sums[assigned_strand] / float(tx.length)
        tx.attrs[DENSITY] = density
        tx.attrs[RECUR] = recur

def annotate_transcripts(known_genes_bed_file, lib_table_file, 
                         transcripts_gtf_file, output_dir):
    # read library table
    logging.debug("Reading library file '%s'" % (lib_table_file))
    libinfos = list(LibInfo.from_file(lib_table_file))
    # make sample -> library dict
    sample_lib_dict = collections.defaultdict(lambda: set())
    for libinfo in libinfos:
        sample_lib_dict[libinfo.sample].add(libinfo.library)
    # compute # of libraries per sample
    libs_per_sample = {}
    for sample,libs in sample_lib_dict.iteritems():
        libs_per_sample[sample] = len(libs)
    # read known transcripts
    logging.debug("Reading known genes file '%s'" % (known_genes_bed_file))
    exon_trees, intron_trees  = build_known_transcript_trees(open(known_genes_bed_file))
    # create directory for classification tree files
    ctree_dir = os.path.join(output_dir, "ctrees")
    if not os.path.exists(ctree_dir):
        logging.debug("Creating directory '%s'" % (ctree_dir))
        os.makedirs(ctree_dir)
    # read transcripts by locus
    logging.debug("Annotating transcripts")
    header = ["chrom", "start", "tx_id", "annotated", "sense_frac", "antisense_frac", 
              "intron_frac", "length", "num_exons", "cov", "fpkm", 
              "recur", "avgdensity"] 
    library_fh_dict = {}
    library_filename_dict = {}
    for locus_transcripts in parse_gtf(open(transcripts_gtf_file)):
        # compute recurrence and expression density
        annotate_recurrence_and_density(locus_transcripts, libs_per_sample)        
        # compute overlap with known genes
        annotate_overlap(locus_transcripts, exon_trees, intron_trees)
        for tx in locus_transcripts:
            # define "annotated" transcripts that have exonic overlap
            # greater than a pre-specified cutoff
            annotated = int(tx.attrs[SENSE_EXON_FRAC] >= MIN_SENSE_EXON_FRAC)
            fields = [tx.chrom,
                      tx.start,
                      tx.id,
                      annotated,
                      '%.3f' % tx.attrs[SENSE_EXON_FRAC],
                      '%.3f' % tx.attrs[ANTISENSE_EXON_FRAC],
                      '%.3f' % tx.attrs[INTRON_FRAC],
                      tx.length,
                      len(tx.exons),
                      tx.cov,
                      tx.fpkm,
                      tx.attrs[RECUR],
                      tx.attrs[DENSITY]]
            # write transcript to file
            if tx.library not in library_fh_dict:
                filename = os.path.join(ctree_dir, tx.library + ".txt")
                library_filename_dict[tx.library] = filename
                fh = open(filename, "w")
                library_fh_dict[tx.library] = fh
                print >>fh, '\t'.join(header)
            fh = library_fh_dict[tx.library]                
            print >>fh, '\t'.join(map(str, fields))
    for fh in library_fh_dict.itervalues():
        fh.close()
    return library_filename_dict

def classify_transcripts(library_filename_dict, result_file, cutoff_file):
    # open output files
    result_fh = open(result_file, "w")
    cutoff_fh = open(cutoff_file, "w")
    for library, filename in library_filename_dict.iteritems():
        # setup classification result files on per-library basis
        prefix = os.path.splitext(filename)[0]
        outfile = prefix + ".classify.txt"
        cutfile = prefix + ".classify.cutoff"
        plotfile = prefix + ".classify.pdf"
        logfile = prefix + ".classify.log"
        logfh = open(logfile, "w")
        # run R script to do classification
        retcode = subprocess.call(["Rscript", R_SCRIPT, filename, outfile, 
                                   cutfile, plotfile], stdout=logfh, 
                                   stderr=logfh)
        logfh.close()
        if retcode != 0:
            logging.error("\t%s - FAILED classification" % (library))
            sys.exit(1)
        # store cutoff
        logging.debug("\t%s finished classification" % (library))
        cutoff = float(open(cutfile).next().strip())
        print >>cutoff_fh, '\t'.join(map(str, [library, cutoff]))
        # append results to master output file
        shutil.copyfileobj(open(outfile), result_fh)
    result_fh.close()
    cutoff_fh.close()

def sort_classification_results(input_file, output_file, tmp_dir):
    # sort classification results
    def sort_by_chrom_start(line):
        fields = line.strip().split('\t', 2)
        return fields[0], int(fields[1])
    batch_sort(input=input_file,
               output=output_file,
               key=sort_by_chrom_start,
               buffer_size=(1 << 21),
               tempdirs=[tmp_dir])

class ClassificationResult(object):
    __slots__ = ('chrom', 'start', 'tx_id', 'annotated', 'sense_frac', 
                 'antisense_frac', 'intron_frac', 'recur', 
                 'avgdensity', 'prob')
    @staticmethod
    def from_line(line):
        fields = line.strip().split('\t')
        c = ClassificationResult()
        c.chrom = fields[0]
        c.start = int(fields[1])
        c.tx_id = fields[2]
        c.annotated = int(fields[3])
        c.sense_frac = float(fields[4])
        c.antisense_frac = float(fields[5])
        c.intron_frac = float(fields[6])
        c.recur = float(fields[7])
        c.avgdensity = float(fields[8])
        c.prob = float(fields[9])
        return c

def summarize_results(result_file, cutoff_file, gtf_file, 
                      expressed_file, background_file):
    # load cutoff file
    library_cutoff_dict = collections.defaultdict(lambda: 1.0)
    for line in open(cutoff_file):
        fields = line.strip().split('\t')
        library_cutoff_dict[fields[0]] = float(fields[1])
    # setup output files
    expressed_fh = open(expressed_file, "w")
    background_fh = open(background_file, "w")
    # keep track of prediction outcomes
    tp, fp, tn, fn = 0, 0, 0, 0
    attr_names = ('annotated', 'sense_frac', 'antisense_frac', 'intron_frac', 'recur', 'avgdensity', 'prob')
    # maintain predictions in a dictionary where key is transcript id and 
    # value is boolean prediction decision for transcript.
    tx_id_decisions = {}
    # maintain heap queues that keeps track of the last transcript position 
    # on each chromosome. prediction decisions only need to be remembered
    # until the parsing goes past the end of the transcript (all exons 
    # accounted for)
    decision_heapqs = collections.defaultdict(lambda: [])
    # maintain result objects in dictionary keyed by transcript id
    tx_id_results = {}
    # maintain heapq that keeps track of transcript position of
    # result objects.  results only need to be remembered until parsing
    # goes past the chrom/start position of the transcript (all transcripts
    # accounted for)
    result_heapqs = collections.defaultdict(lambda: [])
    # read result file and gtf file in sync
    result_fh = open(result_file)
    gtf_fh = open(gtf_file)
    for feature in gtf.GTFFeature.parse(gtf_fh):
        # get transcript id used to lookup expressed/background 
        # prediction decision
        tx_id = feature.attrs["transcript_id"]
        # check top of the decision heapq and pop transcripts when parsing 
        # has gone past the end
        decision_heapq = decision_heapqs[feature.seqid]               
        while (len(decision_heapq) > 0) and (feature.start > decision_heapq[0][0]):
            smallest_end, smallest_tx_id = heapq.heappop(decision_heapq)
            del tx_id_decisions[smallest_tx_id]
        # check top of result heapq and pop transcripts when parsing has gone
        # past the end
        result_heapq = result_heapqs[feature.seqid]
        while (len(result_heapq) > 0) and (feature.start > result_heapq[0][0]):
            result_start, result_tx_id = heapq.heappop(result_heapq)
            del tx_id_results[result_tx_id]
        # parse transcript/exon features differently
        if feature.feature_type == "transcript":
            # push transcript end onto decision heap queue 
            # (decision must stay valid until past the end)
            heapq.heappush(decision_heapq, (feature.end, tx_id))
            # parse results until this tx_id is found (all results 
            # must stay valid until past this chrom/start location)
            while tx_id not in tx_id_results:
                result = ClassificationResult.from_line(result_fh.next())
                assert (result.chrom == feature.seqid)                
                assert (result.start <= feature.start)
                # add to result dictionary
                tx_id_results[result.tx_id] = result
                # add to heapq to remove results that are no longer useful
                heapq.heappush(result_heapqs[result.chrom], (result.start, result.tx_id))
            # lookup classification result
            result = tx_id_results[tx_id]
            # ensure that transcript_id attribute matches result id
            assert tx_id == result.tx_id
            # get library name from gene_id attribute (by convention)
            library = feature.attrs["gene_id"].split(".")[0]
            # lookup cutoff value for classification
            cutoff = library_cutoff_dict[library]                
            # add results as transcript attributes
            for attr_name in attr_names:
                feature.attrs[attr_name] = getattr(result, attr_name)
            # keep track of prediction statistics
            if result.annotated and (result.prob > cutoff):
                tp += 1
            elif result.annotated:
                fn += 1
            elif (result.prob > cutoff):
                fp += 1
            else:
                tn += 1
            # remember decision so that it can be applied to the
            # transcript exons as well
            if result.annotated or (result.prob > cutoff):
                tx_id_decisions[tx_id] = True
                print >>expressed_fh, str(feature)
            else:
                tx_id_decisions[tx_id] = False
                print >>background_fh, str(feature)
        elif feature.feature_type == "exon":
            # get the decision associated with the exon's transcript id
            if tx_id_decisions[tx_id]:
                print >>expressed_fh, str(feature)
            else:
                print >>background_fh, str(feature)
    # cleanup
    result_fh.close()    
    expressed_fh.close()
    background_fh.close()
    # compute statistics
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
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    # parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--output-dir", dest="output_dir", default="transcripts") 
    parser.add_argument("known_genes_file")
    parser.add_argument("lib_table")
    parser.add_argument("gtf_file")
    args = parser.parse_args()

    # setup output dir
    if not os.path.exists(args.output_dir):
        logging.debug("Creating output directory '%s'" % (args.output_dir))
        os.makedirs(args.output_dir)
    # setup output file names
    result_file = os.path.join(args.output_dir, "classification_results.txt")
    cutoff_file = os.path.join(args.output_dir, "classification_cutoffs.txt")
    sorted_result_file = os.path.join(args.output_dir, "classification_results.srt.txt")
    expressed_file = os.path.join(args.output_dir, "expressed.gtf")
    background_file = os.path.join(args.output_dir, "background.gtf")

    logging.info("Annotating transcripts")
    library_filename_dict = annotate_transcripts(args.known_genes_file, args.lib_table, 
                                                 args.gtf_file, args.output_dir)
    logging.info("Classifying transcripts")
    classify_transcripts(library_filename_dict, result_file, cutoff_file)
    logging.info("Sorting classification results")
    sort_classification_results(result_file, sorted_result_file, args.output_dir)
    logging.info("Summarizing results")
    summarize_results(sorted_result_file, cutoff_file, 
                      args.gtf_file, expressed_file, background_file)
    
if __name__ == '__main__': 
    main()