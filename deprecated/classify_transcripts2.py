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
from assemblyline.lib.transcript_parser import parse_gtf, cufflinks_attr_defs
from assemblyline.lib.sampletable import LibInfo
from assemblyline.lib.bx.cluster import ClusterTree
from assemblyline.lib.bx.intersection import Interval, IntervalTree
from assemblyline.lib.batch_sort import batch_sort
from assemblyline.lib import bed
from assemblyline.lib import gtf

# constant transcript attribute names
DENSITY = 'd'
RECUR = 'r'
SENSE_EXON_FRAC = 'sense_frac'
ANTISENSE_EXON_FRAC = 'antisense_frac'
INTRON_FRAC = 'intron_frac'

# constant category names
INTERGENIC = 0
WITHIN_GENES = 1
category_int_to_str = {INTERGENIC: "intergenic", WITHIN_GENES: "within_gene"}
category_str_to_int = dict((v,k) for k,v in category_int_to_str.items())

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

def categorize_transcripts(known_genes_bed_file, transcripts_gtf_file,
                           output_gtf_file_dict):
    # read known transcripts
    logging.debug("Reading known genes file '%s'" % (known_genes_bed_file))
    exon_trees, intron_trees = build_known_transcript_trees(open(known_genes_bed_file))
    # setup output gtf files by category
    fh_dict = {}
    for category, filename in output_gtf_file_dict.iteritems():
        fh_dict[category] = open(filename, "w")
    # read transcripts by locus
    logging.debug("Categorizing transcripts")
    for locus_transcripts in parse_gtf(open(transcripts_gtf_file), cufflinks_attr_defs):
        # compute overlap with known genes
        annotate_overlap(locus_transcripts, exon_trees, intron_trees)
        category_features = {INTERGENIC: [], WITHIN_GENES: []}
        for tx in locus_transcripts:
            # define "annotated" transcripts that have exonic overlap
            # greater than a pre-specified cutoff
            annotated = int(tx.attrs[SENSE_EXON_FRAC] >= MIN_SENSE_EXON_FRAC)
            # categorize unannotated transcript as 
            # 'within_gene' (or intergenic)
            if ((tx.attrs[SENSE_EXON_FRAC] > 0) or
                (tx.attrs[INTRON_FRAC] > 0) or
                (tx.attrs[ANTISENSE_EXON_FRAC] > 0)):
                # within genes category
                category = WITHIN_GENES
            else:
                category = INTERGENIC
            # add to transcript attributes
            tx.attrs["annotated"] = annotated
            tx.attrs["category"] = category_int_to_str[category]
            # convert from Transcript object to GTF
            features = tx.to_gtf_features()
            if not annotated:
                category_features[category].extend(features)
            else:
                category_features[INTERGENIC].extend(features)
                category_features[WITHIN_GENES].extend(features)
        # sort and output category features
        for category, features in category_features.iteritems():
            # sort so that 'transcript' features appear before 'exon'
            features.sort(key=lambda f: f.feature_type, reverse=True)
            features.sort(key=lambda f: f.start)
            # output transcripts to gtf
            for f in features:
                print >>fh_dict[category], str(f) 
    # cleanup
    for fh in fh_dict.itervalues():
        fh.close()

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

def annotate_transcripts(lib_table_file, transcripts_gtf_file, ctree_dir, 
                         lib_info_file):
    # read library table
    logging.debug("Reading library file '%s'" % (lib_table_file))
    libinfos = list(LibInfo.from_file(lib_table_file))
    # make sample -> library dict and
    sample_lib_dict = collections.defaultdict(lambda: set())
    for libinfo in libinfos:
        sample_lib_dict[libinfo.sample].add(libinfo.library)
    # compute # of libraries per sample
    libs_per_sample = {}
    for sample,libs in sample_lib_dict.iteritems():
        libs_per_sample[sample] = len(libs)
    # setup library annotation files by creating a dictionary of 
    # open file handles to write data to
    header = ["chrom", "start", "tx_id", "annotated", "length", 
              "num_exons", "cov", "fpkm", "recur", "avgdensity"]
    filename_dict = {}
    fh_dict = {}
    lib_count_dict = collections.defaultdict(lambda: [0, 0])
    for libinfo in libinfos:
        filename = os.path.join(ctree_dir, libinfo.library + ".txt")
        filename_dict[libinfo.library] = filename 
    # read transcripts by locus
    logging.debug("Annotating transcripts")
    for locus_transcripts in parse_gtf(open(transcripts_gtf_file), cufflinks_attr_defs):
        # compute recurrence and expression density
        annotate_recurrence_and_density(locus_transcripts, libs_per_sample)        
        for tx in locus_transcripts:
            # setup list of annotation fields
            fields = [tx.chrom,
                      tx.start,
                      tx.id,
                      tx.attrs["annotated"],
                      tx.length,
                      len(tx.exons),
                      tx.cov,
                      tx.fpkm,
                      tx.attrs[RECUR],
                      tx.attrs[DENSITY]]
            # lookup file handle and open new file if necessary
            if not tx.library in fh_dict:
                fh_dict[tx.library] = open(filename_dict[tx.library], "w")            
                print >>fh_dict[tx.library], '\t'.join(header)
            # write to file            
            print >>fh_dict[tx.library], '\t'.join(map(str, fields))
            # keep track of # of annotated/unannotated transcripts per file
            annotated = int(tx.attrs["annotated"])
            lib_count_dict[tx.library][annotated] += 1
    # close open file handles
    for fh in fh_dict.itervalues():
        fh.close()
    # write library filename and count information
    fh = open(lib_info_file, "w")
    for library,filename in filename_dict.iteritems():
        counts = lib_count_dict[library]
        print >>fh, '\t'.join(map(str, [library, counts[0], counts[1], filename]))
    fh.close()

def classify_transcripts(lib_info_file, result_file, cutoff_file):
    # read library, category, count information to determine
    # which files need to be classified
    classify_task_list = []
    for line in open(lib_info_file):
        fields = line.strip().split('\t')
        library = fields[0]
        unannotated_count = int(fields[1])
        annotated_count = int(fields[2])
        filename = fields[3]
        total_count = unannotated_count + annotated_count
        if total_count == 0:
            logging.debug("\tskipping empty library %s" % (library))
            continue
        classify_task_list.append((library, filename))
    # open output files
    result_fh = open(result_file, "w")
    cutoff_fh = open(cutoff_file, "w")
    for library, filename in classify_task_list:
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
    __slots__ = ('chrom', 'start', 'tx_id', 'annotated', 'recur', 
                 'avgdensity', 'prob')
    @staticmethod
    def from_line(line):
        fields = line.strip().split('\t')
        c = ClassificationResult()
        c.chrom = fields[0]
        c.start = int(fields[1])
        c.tx_id = fields[2]
        c.annotated = int(fields[3])
        c.recur = float(fields[4])
        c.avgdensity = float(fields[5])
        c.prob = float(fields[6])
        return c

def summarize_results(result_file, cutoff_file, gtf_file, 
                      annotated_file, expressed_file, background_file,
                      stats_file):
    # load cutoff file
    library_cutoff_dict = collections.defaultdict(lambda: 1.0)
    for line in open(cutoff_file):
        fields = line.strip().split('\t')
        library = fields[0]
        cutoff = float(fields[1])
        library_cutoff_dict[library] = float(cutoff)
    # keep track of prediction outcomes
    tp, fp, tn, fn = 0, 0, 0, 0
    attr_names = ('recur', 'avgdensity', 'prob')
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
    # open output files
    annotated_fh = open(annotated_file, "w")
    expressed_fh = open(expressed_file, "w")
    background_fh = open(background_file, "w")
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
            # add cutoff to feature
            feature.attrs["prob_cutoff"] = cutoff
            # keep track of prediction decision and statistics
            # remember decision in tx_id_decision dict so that it can be 
            # applied to the transcript exons as well
            if result.annotated:
                tx_id_decisions[tx_id] = "annotated"
                if result.prob > cutoff:
                    tp += 1
                else:
                    fn += 1
            else:
                if result.prob > cutoff:
                    tx_id_decisions[tx_id] = "expressed"
                    fp += 1
                else:
                    tx_id_decisions[tx_id] = "background"
                    tn += 1
        # lookup the decision associated with the transcript id
        decision = tx_id_decisions[tx_id]
        if decision == "annotated":
            print >>annotated_fh, str(feature)
        elif decision == "expressed":
            print >>expressed_fh, str(feature)
        else:
            assert decision == "background"
            print >>background_fh, str(feature)
    # cleanup
    annotated_fh.close()
    expressed_fh.close()
    background_fh.close()
    result_fh.close()    
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
    fh = open(stats_file, "w")
    print >>fh, "TP: %d" % (tp)
    print >>fh, "FN: %d" % (fn)
    print >>fh, "FP: %d" % (fp)
    print >>fh, "TN: %d" % (tn)
    print >>fh, "Precision: %f" % (prec)
    print >>fh, "Recall: %f" % (rec)
    print >>fh, "Sensitivity: %f" % (sens)
    print >>fh, "Specificity: %f" % (spec)
    print >>fh, "F: %f" % (F)
    print >>fh, "Accuracy: %f" % (acc)
    fh.close()
    logging.info("TP=%d FN=%d FP=%d TN=%d" % (tp, fn, fp, tn))
    logging.info("Precision=%f Recall/Sensitivity=%f" % (prec, rec))
    logging.info("Specificity=%f" % (spec))
    logging.info("F=%f" % F)
    logging.info("Accuracy=%f" % (acc))

def merge_sort_gtf_files(gtf_files, output_file):
    tmp_file = os.path.splitext(output_file)[0] + ".unsorted.gtf"
    outfh = open(tmp_file, "w")
    for filename in gtf_files:
        shutil.copyfileobj(open(filename), outfh)
    outfh.close()
    gtf.sort_gtf(tmp_file, output_file)
    os.remove(tmp_file)

    
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
    # setup output file names for each category
    class CategoryData(object):
        pass
    category_data = {}
    for category in (INTERGENIC, WITHIN_GENES):
        c = CategoryData()
        c.output_dir = os.path.join(args.output_dir, category_int_to_str[category])
        if not os.path.exists(c.output_dir):
            logging.debug("Creating directory: %s" % (c.output_dir))
            os.makedirs(c.output_dir)
        c.gtf_file = os.path.join(c.output_dir, "transcripts.srt.gtf")
        c.ctree_dir = os.path.join(c.output_dir, "ctrees")
        c.lib_info_file = os.path.join(c.output_dir, "library_file_info.txt")
        c.result_file = os.path.join(c.output_dir, "classification_results.txt")
        c.cutoff_file = os.path.join(c.output_dir, "classification_cutoffs.txt")
        c.sorted_result_file = os.path.join(c.output_dir, "classification_results.srt.txt")
        c.annotated_file = os.path.join(c.output_dir, "annotated.gtf")
        c.expressed_file = os.path.join(c.output_dir, "expressed.gtf")
        c.background_file = os.path.join(c.output_dir, "background.gtf")
        c.stats_file = os.path.join(c.output_dir, "stats.txt")
        category_data[category] = c

    # place transcripts into categories
    logging.info("Categorizing transcripts")
    category_gtf_files = {INTERGENIC: category_data[INTERGENIC].gtf_file,
                          WITHIN_GENES: category_data[WITHIN_GENES].gtf_file}
    categorize_transcripts(args.known_genes_file, args.gtf_file, category_gtf_files)

    # merge annotated/expressed transcripts into a single GTF
    merged_expressed_gtf_files = [category_data[INTERGENIC].annotated_file]
    merged_background_gtf_files = []

    # process transcripts in each category         
    for category in (INTERGENIC, WITHIN_GENES):
        logging.info("Classifying category %s" % (category_int_to_str[category]))
        c = category_data[category]
        if not os.path.exists(c.ctree_dir):
            logging.debug("\tcreating directory: %s" % (c.ctree_dir))
            os.makedirs(c.ctree_dir)
        logging.info("\tannotating transcripts gtf_file=%s output_dir=%s" % (c.gtf_file, c.ctree_dir))
        annotate_transcripts(args.lib_table, c.gtf_file, c.ctree_dir, c.lib_info_file) 
        logging.info("\tclassifying transcripts")
        classify_transcripts(c.lib_info_file, c.result_file, c.cutoff_file)
        logging.info("\tsorting classification results")
        sort_classification_results(c.result_file, c.sorted_result_file, c.output_dir)
        logging.info("\tsummarizing results")
        summarize_results(c.sorted_result_file, c.cutoff_file, 
                          c.gtf_file, c.annotated_file, 
                          c.expressed_file, c.background_file,
                          c.stats_file)
        merged_expressed_gtf_files.append(c.expressed_file)
        merged_background_gtf_files.append(c.background_file)

    # merge the GTF files from different categories
    logging.info("Merging classification results")
    expressed_gtf_filename = os.path.join(args.output_dir, "expressed.gtf")
    background_gtf_filename = os.path.join(args.output_dir, "background.gtf")
    merge_sort_gtf_files(merged_expressed_gtf_files, expressed_gtf_filename)
    merge_sort_gtf_files(merged_background_gtf_files, background_gtf_filename)
    
if __name__ == "__main__":
    main()
