'''
Created on Nov 30, 2011

@author: mkiyer
'''
import logging
import argparse
import os
import collections
import numpy as np

from assemblyline.lib.transcript import NO_STRAND, POS_STRAND, NEG_STRAND
from assemblyline.lib.transcript_parser import parse_gtf
from assemblyline.lib.sampletable import LibInfo
from assemblyline.lib.bx.cluster import ClusterTree
from assemblyline.lib.bx.intersection import Interval, IntervalTree
from assemblyline.lib import bed

# constant transcript attribute names
DENSITY = 'd'
RECUR = 'r'
SENSE_EXON_FRAC = 'sf'
ANTISENSE_EXON_FRAC = 'asf'
INTRON_FRAC = 'if'

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
    locus_end = transcripts[-1].end
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

    # compute recurrence/density of individual transcripts using the arrays
    for tx in transcripts:
        # sum mass/recurrence on all strands over the exonic intervals
        # of the transcript
        recur_sums = np.zeros(3, dtype=int)
        mass_sums = np.zeros(3, dtype=float)        
        for e in tx.exons:           
            recur_sums += np.sum(recur_arr[:,(e.start - locus_start):(e.end - locus_start)], axis=1)
            mass_sums += np.sum(mass_arr[:,(e.start - locus_start):(e.end - locus_start)], axis=1)
        # assign most likely strand to unstranded transcripts
        total_stranded_mass = mass_sums[POS_STRAND] + mass_sums[NEG_STRAND]
        if tx.strand == NO_STRAND:
            assigned_strand = np.argmax(mass_sums)
        else:
            assigned_strand = tx.strand
        # calculate strand fraction
        if total_stranded_mass == 0:
            # all of the mass is 'unstranded'
            strand_frac = 1.0
        else:
            strand_frac = mass_sums[assigned_strand] / total_stranded_mass
        # calculate recurrence and density (mass/length)
        density = (mass_sums[assigned_strand] + strand_frac * mass_sums[NO_STRAND]) / float(tx.length)
        recur = (recur_sums[assigned_strand] + strand_frac * recur_sums[NO_STRAND]) / float(tx.length)
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
    header = ["tx_id", "known", "within_gene", "length", "num_exons",
              "cov", "fpkm", "recur", "avgdensity"] 
    library_fh_dict = {}
    for locus_transcripts in parse_gtf(open(transcripts_gtf_file)):
        # compute recurrence and expression density
        annotate_recurrence_and_density(locus_transcripts, libs_per_sample)        
        # compute overlap with known genes
        annotate_overlap(locus_transcripts, exon_trees, intron_trees)
        for tx in locus_transcripts:
            if tx.attrs[SENSE_EXON_FRAC] > 0:
                # overlaps known gene
                known = 1
                within_gene = 1
            elif ((tx.attrs[INTRON_FRAC] > 0) or
                  (tx.attrs[ANTISENSE_EXON_FRAC] > 0)):
                # within gene boundaries
                known = 0
                within_gene = 1
            else:
                # intergenic
                known = 0
                within_gene = 0
            fields = [tx.id,
                      known,
                      within_gene,
                      tx.length,
                      len(tx.exons),
                      tx.cov,
                      tx.fpkm,
                      tx.attrs[RECUR],
                      tx.attrs[DENSITY]]
            if tx.library not in library_fh_dict:
                filename = os.path.join(ctree_dir, tx.library + ".txt")
                fh = open(filename, "w")
                library_fh_dict[tx.library] = fh
                print >>fh, '\t'.join(header)
            fh = library_fh_dict[tx.library]                
            print >>fh, '\t'.join(map(str, fields))

def classify_transcripts(output_dir):
    # setup output files
    expressed_file = os.path.join(output_dir, "expressed.gtf")
    background_file = os.path.join(output_dir, "background.gtf")
    expressed_fh = open(expressed_file, "w")
    background_fh = open(background_file, "w")
    # cleanup    
    expressed_fh.close()
    background_fh.close()

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
    annotate_transcripts(args.known_genes_file, args.lib_table, 
                         args.gtf_file, args.output_dir)
    #classify_transcripts(args.known_genes_file, args.lib_table, args.gtf_file, args.output_dir)
    
if __name__ == '__main__': 
    main()