'''
Created on Jun 5, 2010
Modified on Mar 29, 2011

@author: mkiyer
'''
import argparse
import collections
import logging
import os

from bx.intervals.intersection import Interval, IntervalTree
from bx.intervals.cluster import ClusterTree
import lib.bed as bed

def parse_tabular(line,
                  chrom_col=0,
                  start_col=1,
                  end_col=2,
                  strand_col=3):
    if line == None:
        return None
    line = line.strip()
    if line.startswith('#'):
        return None
    fields = line.split('\t')
    return Interval(int(fields[start_col]),
                    int(fields[end_col]),
                    chrom=fields[chrom_col],
                    strand=fields[strand_col])

def build_exon_trees(bedfile, padding):
    trees = collections.defaultdict(lambda: IntervalTree())
    for g in bed.BEDFeature.parse(open(bedfile)):
        for start, end in g.exons:
            trees[g.chrom].insert_interval(Interval(start-padding, end+padding, strand=g.strand))
    return trees

def build_interval_trees(bedfile, padding):
    gene_trees = collections.defaultdict(lambda: IntervalTree())
    for line in open(bedfile):
        fields = line.split('\t')
        chrom = fields[0]
        start = int(fields[1])
        end = int(fields[2])
        strand = fields[5]
        gene_trees[chrom].insert_interval(Interval(start-padding, end+padding, strand=strand))    
    return gene_trees

def cmp_strand(a, b):
    if a == b:
        return True
    if (a == '.') or (b == '.'):
        return True
    return False

def intersect_bed12(bedfile, format_args, trees, exons=True, stranded=True):
    for line in open(bedfile):
        g = bed.BEDFeature.from_string(line)
        result = False
        if exons:
            for start, end in g.exons:
                for hit in trees[g.chrom].find(start, end):
                    if stranded and (not cmp_strand(g.strand, hit.strand)):
                        continue
                    result = True
                    break
                if result:
                    break
        else:
            for hit in trees[g.chrom].find(g.tx_start, g.tx_end):
                if stranded and (not cmp_strand(g.strand, hit.strand)):
                    continue
                result = True
                break            
        yield result, line

def intersect_tabular(infile, format_args, trees, exons, stranded):
    chrom_col = format_args['chrom_col']
    start_col = format_args['start_col']
    end_col = format_args['end_col']
    strand_col = format_args['strand_col']
    for line in open(infile):
        interval = parse_tabular(line, chrom_col, start_col, end_col, strand_col)
        result = False
        for hit in trees[interval.chrom].find(interval.start, interval.end):
            if stranded and (not cmp_strand(interval.strand, hit.strand)):
                continue
            result = True
            break            
        yield result, line


def intersect(infile, format, format_args, trees, exons=True, stranded=True):
    if format == "bed12":
        func = intersect_bed12
    elif format == "tabular":
        func = intersect_tabular
    for result, line in func(infile, format_args, trees, exons, stranded):
        yield result, line

def do_intersection(input_file, format, format_args, hitfile, missfile, trees, exons, stranded):
    hitfileh = open(hitfile, "w")
    missfileh = open(missfile, "w")
    # intersect intervals
    for hit, line in intersect(input_file, format, format_args, trees, exons, stranded):
        if hit:
            hitfileh.write(line)
        else:
            missfileh.write(line)
    hitfileh.close()
    missfileh.close()

def categorize_transcripts(bedfile, format, format_args, basename, output_dir,
                           pseudo_bed_file, protein_bed_file, 
                           annotated_bed_file):
    # build interval trees
    logging.debug("Loading pseudogenes...")
    pseudo_trees = build_interval_trees(pseudo_bed_file, 0)
    logging.debug("Loading protein coding gene exons...")
    protein_exon_trees = build_exon_trees(protein_bed_file, 0)
    logging.debug("Loading protein coding genes...")
    protein_gene_trees = build_interval_trees(protein_bed_file, 0)
    logging.debug("Loading annotated exons...")
    annotated_exon_trees = build_exon_trees(annotated_bed_file, 0)
    #logging.debug("Loading annotated genes...")
    #annotated_gene_trees = build_interval_trees(annotated_bed_file, 0)
    logging.debug("Categorizing transcripts...")
    # pseudo genes
    input_file = bedfile
    pseudo_file = os.path.join(output_dir, basename + '.pseudo.bed') 
    notpseudo_file = os.path.join(output_dir, basename + '.notpseudo.bed')
    # intersect with pseudogenes
    hit_file, miss_file = pseudo_file, notpseudo_file
    do_intersection(input_file, format, format_args, hit_file, miss_file, pseudo_trees, exons=True, stranded=False)

    # annotated genes (not pseudogenes)
    input_file = notpseudo_file
    prefix = os.path.splitext(input_file)[0]
    annot_file = os.path.join(output_dir, prefix + '.annotated.bed')
    notannot_file = os.path.join(output_dir, prefix + '.notannotated.bed')
    # intersect non-pseudogenes with annotated genes
    hit_file, miss_file = annot_file, notannot_file
    do_intersection(input_file, format, format_args, hit_file, miss_file, annotated_exon_trees, exons=True, stranded=True)

    # completely unannotated (including strand)
    input_file = notpseudo_file
    prefix = os.path.splitext(input_file)[0]
    nostrand_annot_file = os.path.join(output_dir, prefix + '.nostrand_annotated.bed')
    nostrand_notannot_file = os.path.join(output_dir, prefix + '.nostrand_notannotated.bed')
    # intersect non-pseudogenes with annotated genes
    hit_file, miss_file = nostrand_annot_file, nostrand_notannot_file
    do_intersection(input_file, format, format_args, hit_file, miss_file, annotated_exon_trees, exons=True, stranded=False)

    # annotated proteins
    input_file = annot_file
    prefix = os.path.splitext(input_file)[0]
    annot_protein_file = os.path.join(output_dir, prefix + '.protein.bed') 
    annot_notprotein_file = os.path.join(output_dir, prefix + '.notprotein.bed')
    # intersect annotated genes with protein coding genes
    hit_file, miss_file = annot_protein_file, annot_notprotein_file
    do_intersection(input_file, format, format_args, hit_file, miss_file, protein_exon_trees, exons=True, stranded=True)

    # annotated antisense non-coding RNAs
    input_file = annot_notprotein_file
    prefix = os.path.splitext(input_file)[0]
    annot_notprotein_antisense_file = os.path.join(output_dir, prefix + '.antisense.bed') 
    annot_notprotein_notantisense_file = os.path.join(output_dir, prefix + '.notantisense.bed') 
    # intersect
    hit_file, miss_file = annot_notprotein_antisense_file, annot_notprotein_notantisense_file
    do_intersection(input_file, format, format_args, hit_file, miss_file, protein_exon_trees, exons=True, stranded=False)

    # annotated non-coding RNAs
    input_file = annot_notprotein_notantisense_file
    prefix = os.path.splitext(input_file)[0]
    annot_notprotein_intronic_file = os.path.join(output_dir, prefix + '.intronic.bed') 
    annot_notprotein_intergenic_file = os.path.join(output_dir, prefix + '.intergenic.bed')
    # intersect annotated non-codings to see whether the are intronic/intergenic
    hit_file, miss_file = annot_notprotein_intronic_file, annot_notprotein_intergenic_file
    do_intersection(input_file, format, format_args, hit_file, miss_file, protein_gene_trees, exons=False, stranded=False)

    # unannotated antisense transcripts that overlap proteins
    input_file = notannot_file
    prefix = os.path.splitext(input_file)[0]
    notannot_antisense_file = os.path.join(output_dir, prefix + '.antisense.bed') 
    notannot_notantisense_file = os.path.join(output_dir, prefix + '.notantisense.bed') 
    # intersect
    hit_file, miss_file = notannot_antisense_file, notannot_notantisense_file
    do_intersection(input_file, format, format_args, hit_file, miss_file, protein_exon_trees, exons=True, stranded=False)
    
    # unannotated intronic/intergenic transcripts
    input_file = notannot_notantisense_file
    prefix = os.path.splitext(input_file)[0]
    notannot_intronic_file = os.path.join(output_dir, prefix + '.intronic.bed') 
    notannot_intergenic_file = os.path.join(output_dir, prefix + '.intergenic.bed') 
    # intersect
    hit_file, miss_file = notannot_intronic_file, notannot_intergenic_file
    do_intersection(input_file, format, format_args, hit_file, miss_file, protein_gene_trees, exons=False, stranded=False)
    
    # overall number of protein vs. non-coding
    input_file = notpseudo_file
    prefix = os.path.splitext(input_file)[0]    
    protein_file = os.path.join(output_dir, prefix + '.protein.bed') 
    notprotein_file = os.path.join(output_dir, prefix + '.notprotein.bed') 
    # intersect
    hit_file, miss_file = protein_file, notprotein_file
    do_intersection(input_file, format, format_args, hit_file, miss_file, protein_exon_trees, exons=True, stranded=True)

    # overall antisense noncoding
    input_file = notprotein_file
    prefix = os.path.splitext(input_file)[0]
    notprotein_antisense_file = os.path.join(output_dir, prefix + '.antisense.bed') 
    notprotein_notantisense_file = os.path.join(output_dir, prefix + '.notantisense.bed') 
    # intersect
    hit_file, miss_file = notprotein_antisense_file, notprotein_notantisense_file
    do_intersection(input_file, format, format_args, hit_file, miss_file, protein_exon_trees, exons=True, stranded=False)

    # overall intronic/intergenic noncoding
    input_file = notprotein_notantisense_file
    prefix = os.path.splitext(input_file)[0]
    notprotein_intronic_file = os.path.join(output_dir, prefix + '.intronic.bed') 
    notprotein_intergenic_file = os.path.join(output_dir, prefix + '.intergenic.bed') 
    # intersect
    hit_file, miss_file = notprotein_intronic_file, notprotein_intergenic_file
    do_intersection(input_file, format, format_args, hit_file, miss_file, protein_gene_trees, exons=False, stranded=False)

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--format', dest="format", default="bed12")
    parser.add_argument("--start", dest="start_col", type=int, default=1)
    parser.add_argument("--end", dest="end_col", type=int, default=2)
    parser.add_argument("--chrom", dest="chrom_col", type=int, default=0)
    parser.add_argument("--strand", dest="strand_col", type=int, default=4)    
    parser.add_argument("--pseudo", dest="pseudo_bed_file", required=True)
    parser.add_argument("--protein", dest="protein_bed_file", required=True)
    parser.add_argument("--annotated", dest="annotated_bed_file", required=True)
    parser.add_argument("bedfile")
    parser.add_argument("basename")
    parser.add_argument("output_dir")
    options = parser.parse_args()    
    # check that intersection files exist
    if not os.path.exists(options.pseudo_bed_file):
        parser.error("pseudogene BED file not found")
    if not os.path.exists(options.protein_bed_file):
        parser.error("protein BED file not found")
    if not os.path.exists(options.annotated_bed_file):
        parser.error("annotated BED file not found")
    # run the categorization
    if options.format == "bed12":
        categorize_transcripts(options.bedfile, options.format, {}, options.basename, options.output_dir)
    elif options.format == "tabular":
        format_args = {'chrom_col': options.chrom_col,
                       'start_col': options.start_col,
                       'end_col': options.end_col,
                       'strand_col': options.strand_col}
        categorize_transcripts(options.bedfile, options.format, format_args, options.basename, options.output_dir)

if __name__ == '__main__':
    main()