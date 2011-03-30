'''
Created on Mar 13, 2011

@author: mkiyer
'''
import logging
import argparse
import collections
import os
import sys
import subprocess

from lib.bx.intersection import Interval, IntervalTree
from lib.bx.cluster import ClusterTree
from lib.transcript_parser import parse_gtf
from lib.base import POS_STRAND, NEG_STRAND, NO_STRAND
import lib.bed as bed
import lib.gtf as gtf

# R script for classifying transcripts
R_SCRIPT = os.path.join(os.path.dirname(__file__), "lib", "classify_transcripts.R")

# defaults
DEFAULT_INTRON_OVERLAP_THRESHOLD = 0.9

# transcript attributes
OVERLAP_WEIGHT = 'owt'
OVERLAP_LENGTH = 'olen'
OVERLAP_NUM_TXS = 'ontx'
OVERLAP_RECUR = 'orecur'
OVERLAP_DENSITY = 'odensity'
SENSE_EXON_FRAC = 'sense_exon'
ANTISENSE_EXON_FRAC = 'asense_exon'
INTRON_FRAC = 'intron'


def build_known_transcript_trees(bed_iter):
    exon_trees = {"+": collections.defaultdict(lambda: IntervalTree()),
                  "-": collections.defaultdict(lambda: IntervalTree())}
    intron_trees = collections.defaultdict(lambda: IntervalTree())
    #intergenic_trees = collections.defaultdict(lambda: IntervalTree())
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
        # insert regions between loci into intergenic trees
        #for i1,i2 in zip(loci[:-1], loci[1:]):
        #    intergenic_trees[chrom].insert_interval(Interval(i1[1], i2[0]))
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

def get_overlapping_bases(tx, start, end):
    count = 0
    for e in tx.exons:
        if (e.start >= end) or (e.end <= start):
            continue
        sorted_coords = sorted((e.start, start, e.end, end))
        count += (sorted_coords[2] - sorted_coords[1])
    return count

class StrandData(object):
    def __init__(self):
        self.weight = 0
        self.length = 0
        self.density = 0
        self.recur = set()
        self.num_txs = 0

def annotate_recurrence_and_density(transcripts):
    cluster_tree = ClusterTree(0, 1)
    for f_id, f in enumerate(transcripts):        
        for e in f.exons:
            cluster_tree.insert(e.start, e.end, f_id)
    for start, end, f_ids in cluster_tree.getregions():
        strand_data = [StrandData(), StrandData(), StrandData()]
        # keep track of data on each strand        
        for f_id in f_ids:
            t = transcripts[f_id]
            overlap = get_overlapping_bases(t, start, end)
            score = t.score * (overlap / float(t.length))
            strand_data[t.strand].weight += score
            strand_data[t.strand].length += overlap
            strand_data[t.strand].recur.add(t.label)
            strand_data[t.strand].num_txs += 1
        # compute density
        for strand,data in enumerate(strand_data):
            if data.length == 0:
                data.density = 0
            else:
                data.density = data.weight / float(data.length)
            data.recur = len(data.recur)
        total_density = (strand_data[POS_STRAND].density +
                         strand_data[NEG_STRAND].density)
        if total_density == 0:
            # if there is no "stranded" coverage at this node, 
            # then assign coverage to the positive strand by convention
            pos_frac = 1.0
        else:
            # proportionally assign unstranded coverage based on amount of
            # plus and minus strand coverage
            pos_frac = strand_data[POS_STRAND].density / float(total_density)
        # allocate recurrence/density from unstranded transcripts        
        strand_fracs = (pos_frac, 1.0 - pos_frac)
        for s in (POS_STRAND, NEG_STRAND):
            d = strand_data[s]
            d.weight += strand_fracs[s] * strand_data[NO_STRAND].weight
            d.length += strand_fracs[s] * strand_data[NO_STRAND].length
            d.num_txs += strand_fracs[s] * strand_data[NO_STRAND].num_txs
            d.recur += strand_fracs[s] * strand_data[NO_STRAND].recur
            if d.length == 0:
                d.density = 0
            else:
                d.density = d.weight / float(d.length)
        # assign 'unstranded' transcript stats to the 'best' strand
        d = strand_data[NO_STRAND]
        d.weight = max(strand_data[x].weight for x in (POS_STRAND, NEG_STRAND))               
        d.length = max(strand_data[x].length for x in (POS_STRAND, NEG_STRAND))               
        d.density = max(strand_data[x].density for x in (POS_STRAND, NEG_STRAND))               
        d.num_txs = max(strand_data[x].num_txs for x in (POS_STRAND, NEG_STRAND))               
        d.recur = max(strand_data[x].recur for x in (POS_STRAND, NEG_STRAND))               
        # annotate transcripts with statistics
        for f_id in f_ids:
            tx = transcripts[f_id]            
            if OVERLAP_WEIGHT not in tx.attrs:            
                tx.attrs[OVERLAP_WEIGHT] = 0
                tx.attrs[OVERLAP_LENGTH] = 0
                tx.attrs[OVERLAP_NUM_TXS] = 0
                tx.attrs[OVERLAP_RECUR] = 0
                tx.attrs[OVERLAP_DENSITY] = 0
            d = strand_data[tx.strand]
            tx.attrs[OVERLAP_WEIGHT] += d.weight
            tx.attrs[OVERLAP_LENGTH] += d.length
            tx.attrs[OVERLAP_NUM_TXS] += d.num_txs
            tx.attrs[OVERLAP_RECUR] += d.recur
    # now that all exons are accounted for, compute density 
    for tx in transcripts:        
        tx.attrs[OVERLAP_DENSITY] += tx.attrs[OVERLAP_WEIGHT] / float(tx.attrs[OVERLAP_LENGTH])

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
    return transcripts

def collect_stats(known_genes_file, gtf_file, output_dir,
                  outfh=None, bgfh=None, 
                  intron_overlap_threshold=DEFAULT_INTRON_OVERLAP_THRESHOLD):
    if outfh is None:
        outfh = sys.stdout
    # read known transcripts
    logging.debug("Reading known genes")
    exon_trees, intron_trees  = build_known_transcript_trees(open(known_genes_file))
    # annotate transcripts with properties
    logging.debug("Annotating transcripts")
    header = ["tx_id", "sense_overlap", "antisense_overlap", "intronic_overlap", 
              "length", "num_exons", "mydensity", "recur", "avgdensity"]
    fhdict = {}
    for locus_transcripts in parse_gtf(open(gtf_file)):
        annotate_recurrence_and_density(locus_transcripts)
        annotate_overlap(locus_transcripts, exon_trees, intron_trees)
        for tx in locus_transcripts:            
            fields = [tx.id,
                      tx.attrs[SENSE_EXON_FRAC], 
                      tx.attrs[ANTISENSE_EXON_FRAC],
                      1 if tx.attrs[INTRON_FRAC] >= intron_overlap_threshold else 0,
                      tx.length,
                      len(tx.exons),
                      tx.score / float(tx.length),
                      tx.attrs[OVERLAP_RECUR],
                      tx.attrs[OVERLAP_DENSITY]]
            if tx.label not in fhdict:
                filename = os.path.join(output_dir, tx.label)
                fh = open(filename, "w")
                fhdict[tx.label] = fh
                print >>fh, '\t'.join(header)
            fh = fhdict[tx.label]                
            print >>fh, '\t'.join(map(str, fields))
    # classify transcripts
    logging.debug("Classifying transcripts")
    outdict = {}
    for label, fh in fhdict.iteritems():
        filename = os.path.abspath(fh.name)
        fh.close()
        outfile = os.path.join(output_dir, filename + ".classify.txt")
        cutfile = os.path.join(output_dir, filename + ".classify.cutoff")
        epsfile = os.path.join(output_dir, filename + ".classify.eps")
        logfile = os.path.join(output_dir, filename + ".classify.log")
        logfh = open(logfile, "w")
        retcode = subprocess.call(["Rscript", R_SCRIPT, filename, outfile, 
                                   cutfile, epsfile], stdout=logfh, 
                                   stderr=logfh)
        logfh.close()
        if retcode != 0:
            logging.error("\t%s - FAILED classification" % (label))
        else:
            logging.debug("\t%s" % (label))
            outdict[label] = (outfile, cutfile)
    # write output
    logging.debug("Merging results and writing output")
    # read all of the transcripts
    tx_id_prob_map = {}
    label_cutoff_map = {}
    for label,files in outdict.iteritems():
        outfile, cutfile = files
        for line in open(outfile):
            file_tx_id, prob = line.strip().split('\t')
            prob = float(prob)
            tx_id_prob_map[file_tx_id] = prob
        # get cutoff value
        cutoff = float(open(cutfile).next().strip())
        label_cutoff_map[label] = cutoff
    # write probabilities to GTF
    for feature in gtf.GTFFeature.parse(open(gtf_file)):
        tx_id = feature.attrs["transcript_id"]
        # TODO: need different way to get label
        label = feature.attrs["gene_id"].split(".")[0]
        if tx_id not in tx_id_prob_map:
            logging.warning("\t%s - MISSING data" % (tx_id))
            prob = 0.0
            cutoff = 0.5
        else:
            prob = tx_id_prob_map[tx_id]
            cutoff = label_cutoff_map[label]
        feature.attrs["prob"] = prob
        if prob >= cutoff:
            print >>outfh, str(feature)
        elif bgfh is not None:
            print >>bgfh, str(feature)
    for fh in fhdict.itervalues():
        fh.close()
    

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("--intron-overlap", type=float, 
                        dest="intron_overlap", default=0.90)
    parser.add_argument("--bg", "--background", dest="background_file", 
                        default=None) 
    parser.add_argument("-o", "--expressed", dest="output_file", default=None) 
    parser.add_argument("known_genes_file")
    parser.add_argument("gtf_file")
    parser.add_argument("output_dir")
    options = parser.parse_args()
    # process command line arguments
    output_dir = options.output_dir
    if os.path.exists(output_dir):
        parser.error("output directory already exists.. aborting")        
    logging.debug("Creating output directory")
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)
    if options.output_file is None:
        outfh = None
    else:
        outfh = open(options.output_file, "w")
    if options.background_file is None:
        bgfh = None
    else:
        bgfh = open(options.background_file, "w")    
    collect_stats(options.known_genes_file, options.gtf_file, output_dir,
                  outfh, bgfh, intron_overlap_threshold=options.intron_overlap)
    logging.debug("Finished")

if __name__ == '__main__': main()