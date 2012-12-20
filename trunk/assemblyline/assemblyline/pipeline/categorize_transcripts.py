'''
Created on Jan 25, 2012

@author: mkiyer

AssemblyLine: transcriptome meta-assembly from RNA-Seq

Copyright (C) 2012 Matthew Iyer

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
'''
import logging
import argparse
import collections
import operator

from assemblyline.lib.bx.intersection import Interval, IntervalTree
from assemblyline.lib.bx.cluster import ClusterTree

from assemblyline.lib import gtf
from assemblyline.lib.base import GTFAttr
from assemblyline.lib.transcript import strand_int_to_str
from assemblyline.lib.transcript_parser import parse_gtf, cufflinks_attr_defs

MAX_DIST = 1000000
CATEGORY_PROTEIN = "protein"
CATEGORY_NCRNA = "ncrna"
CATEGORY_ANTISENSE = "antisense"
CATEGORY_INTRONIC = "intronic"
CATEGORY_INTERGENIC = "intergenic"

class Gene(object):
    def __init__(self):
        self.gene_id = None
        self.chrom = None
        self.gene_start = None
        self.gene_end = None
        self.gene_names = set()
        self.exons = set()
        self.strand = "."
        self.is_coding = False
        self.annotation_sources = set()

def cmp_strand(a,b):
    if a == "." or b == ".":
        return True
    return a == b

def get_locus_genes(features):
    gene_map = collections.defaultdict(lambda: Gene())
    for f in features:
        # get gene by id
        gene_id = f.attrs["gene_id"]
        if gene_id not in gene_map:
            g = Gene()
            g.gene_id = gene_id
            g.strand = f.strand
            g.gene_start = f.start
            g.gene_end = f.end
            gene_map[gene_id] = g
        else:
            g = gene_map[gene_id]
        # update gene
        g.gene_start = min(g.gene_start, f.start)
        g.gene_end = max(g.gene_end, f.end)
        if f.feature_type == "exon":
            g.exons.add((f.start, f.end))
        elif f.feature_type == "CDS":
            g.is_coding = True
        if "gene_name" in f.attrs:                    
            g.gene_names.add(f.attrs["gene_name"])
        g.annotation_sources.add(f.source)
    for g in gene_map.itervalues():
        # cluster gene exons
        cluster_tree = ClusterTree(0,1)
        for start,end in g.exons:
            cluster_tree.insert(start, end, 1)
        # update exons
        exon_clusters = []
        for start,end,indexes in cluster_tree.getregions():
            exon_clusters.append((start,end))
        g.exons = exon_clusters
        del cluster_tree
    return sorted(gene_map.values(), key=operator.attrgetter('gene_start'))
    
def read_reference_gtf(ref_gtf_file):
    gene_map = {}
    logging.info("Reading reference GTF file")
    for f in gtf.GTFFeature.parse(open(ref_gtf_file)):
        # get gene by id
        gene_id = f.attrs["gene_id"]
        if gene_id not in gene_map:
            g = Gene()
            g.gene_id = gene_id
            g.chrom = f.seqid 
            g.strand = f.strand
            g.gene_start = f.start
            g.gene_end = f.end
            gene_map[gene_id] = g
        else:
            g = gene_map[gene_id]
        # update gene
        g.gene_start = min(g.gene_start, f.start)
        g.gene_end = max(g.gene_end, f.end)
        if f.feature_type == "exon":
            g.exons.add((f.start, f.end))
        elif f.feature_type == "CDS":
            g.is_coding = True
        if "gene_name" in f.attrs:                    
            g.gene_names.add(f.attrs["gene_name"])
        g.annotation_sources.add(f.source)
    logging.info("Sorting genes")
    genes = sorted(gene_map.values(), key=operator.attrgetter('chrom', 'gene_start'))
    del gene_map
    # cluster loci
    logging.info("Building interval index")
    locus_cluster_trees = collections.defaultdict(lambda: ClusterTree(0,1))
    locus_trees = collections.defaultdict(lambda: IntervalTree())
    for i,g in enumerate(genes):
        locus_cluster_trees[g.chrom].insert(g.gene_start, g.gene_end, i)
    for chrom, cluster_tree in locus_cluster_trees.iteritems(): 
        for locus_start,locus_end,indexes in cluster_tree.getregions():
            # cluster gene exons and add to interval tree
            exon_tree = IntervalTree()
            for i in indexes:
                g = genes[i]
                cluster_tree = ClusterTree(0,1)
                for start,end in g.exons:
                    cluster_tree.insert(start, end, 1)
                # update exons
                exon_clusters = []
                for start,end,indexes in cluster_tree.getregions():
                    exon_clusters.append((start,end))
                g.exons = exon_clusters
                del cluster_tree
                for start,end in g.exons:
                    exon_tree.insert_interval(Interval(start, end, value=g))
            # add to locus interval tree
            locus_trees[chrom].insert_interval(Interval(locus_start, locus_end, value=exon_tree))
    logging.info("Done indexing reference GTF file")
    return locus_trees

def cluster_isoforms(transcripts):
    # cluster exons
    cluster_tree = ClusterTree(0,1)
    for t in transcripts:
        for e in t.exons:
            cluster_tree.insert(e.start, e.end, 1)
    exons = []
    for start,end,indexes in cluster_tree.getregions():
        exons.append((start,end))
    del cluster_tree
    return exons

def get_nearest_genes(chrom, start, end, locus_trees):
    # look left and right
    hits = locus_trees[chrom].before(start, num_intervals=1, max_dist=MAX_DIST)        
    hits.extend(locus_trees[chrom].after(end, num_intervals=1, max_dist=MAX_DIST))
    nearest_locus_hit = None
    nearest_dist = MAX_DIST
    for hit in hits:
        dist = min(abs(start - hit.end), abs(hit.start - end))
        if dist < nearest_dist:
            nearest_dist = dist
            nearest_locus_hit = hit
    nearest_genes = []
    if nearest_locus_hit is None:
        nearest_dist = -1
    else:        
        # find dist to nearest genes
        nearest_dist = MAX_DIST
        gene_tree = nearest_locus_hit.value
        gene_hits = gene_tree.find(nearest_locus_hit.start, nearest_locus_hit.end)
        for gene_hit in gene_hits:
            g = gene_hit.value
            dist = min(abs(start - hit.end), abs(hit.start - end))
            if dist < nearest_dist:
                nearest_dist = dist
        # get nearest gene(s) at distance
        for gene_hit in gene_hits:
            g = gene_hit.value
            dist = min(abs(start - hit.end), abs(hit.start - end))
            if dist == nearest_dist:
                nearest_genes.append(g)
    return nearest_genes, nearest_dist

def categorize_gene_transcripts(transcripts, locus_trees):
    # get single set of exons for all isoforms in gene
    gene_exons = cluster_isoforms(transcripts)
    gene_start = min(e[0] for e in gene_exons)
    gene_end = max(e[1] for e in gene_exons) 
    gene_chrom = transcripts[0].chrom
    gene_strand = strand_int_to_str(transcripts[0].strand)
    # determine whether gene overlaps known loci
    # intersect transcript with reference loci
    locus_hits = locus_trees[gene_chrom].find(gene_start, gene_end)
    if len(locus_hits) == 0:
        # this is a completely unannotated transcript
        category = CATEGORY_INTERGENIC
        nearest_genes, nearest_dist = get_nearest_genes(gene_chrom, gene_start, gene_end, locus_trees)
    else:
        # this transcript overlaps at least one known locus, so
        # categorize as sense/antisense, coding/noncoding, exon/intron
        protein_genes = {}
        ncrna_genes = {}
        antisense_genes = {}
        for locus_hit in locus_hits:
            gene_tree = locus_hit.value
            for exon in gene_exons:
                for gene_hit in gene_tree.find(exon[0], exon[1]):
                    g = gene_hit.value
                    if cmp_strand(g.strand, gene_strand):
                        if g.is_coding:
                            protein_genes[g.gene_id] = g
                        else:
                            ncrna_genes[g.gene_id] = g
                    else:
                        antisense_genes[g.gene_id] = g
        protein_genes = protein_genes.values()
        ncrna_genes = ncrna_genes.values()
        antisense_genes = antisense_genes.values()
        nearest_dist = 0
        if len(protein_genes) > 0:
            category = CATEGORY_PROTEIN
            nearest_genes = protein_genes
        elif len(ncrna_genes) > 0:
            category = CATEGORY_NCRNA
            nearest_genes = ncrna_genes
        elif len(antisense_genes) > 0:
            category = CATEGORY_ANTISENSE
            nearest_genes = antisense_genes
        else:
            category = CATEGORY_INTRONIC
            nearest_genes = []
            for locus_hit in locus_hits:
                gene_tree = locus_hit.value
                gene_hits = gene_tree.find(locus_hit.start, locus_hit.end)
                for gene_hit in gene_hits:
                    nearest_genes.append(gene_hit.value)
    # using 'nearest genes' list get gene names and annotation sources    
    if len(nearest_genes) == 0:
        gene_ids = "NA"
        gene_names = "NA"
        annotation_sources = "NA"
        nearest_dist = -1
    else:
        gene_ids = set()
        gene_names = set()
        annotation_sources = set()
        for g in nearest_genes:
            gene_ids.add(g.gene_id)
            gene_names.update(g.gene_names)
            annotation_sources.update(g.annotation_sources)
        gene_ids = ",".join(sorted(gene_ids))
        gene_names = ",".join(sorted(gene_names))
        annotation_sources = ",".join(sorted(annotation_sources))
    # add attributes to original transcripts
    for t in transcripts:
        t.attrs["category"] = category
        t.attrs["nearest_gene_ids"] = gene_ids
        t.attrs["nearest_gene_names"] = gene_names
        t.attrs["nearest_dist"] = nearest_dist
        t.attrs["annotation_sources"] = annotation_sources

def categorize_transcripts(ref_gtf_file, gtf_file):
    locus_trees = read_reference_gtf(ref_gtf_file)
    logging.debug("Categorizing GTF file")
    for locus_transcripts in parse_gtf(open(gtf_file), attr_defs=cufflinks_attr_defs):
        # group transcripts by gene id
        gene_transcript_map = collections.defaultdict(lambda: [])
        for t in locus_transcripts:
            gene_transcript_map[t.attrs[GTFAttr.GENE_ID]].append(t)
        # categorize genes
        for gene_transcripts in gene_transcript_map.itervalues():
            categorize_gene_transcripts(gene_transcripts, locus_trees)
            # output transcript
            for t in gene_transcripts:
                for f in t.to_gtf_features():
                    print str(f)

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("ref_gtf_file")
    parser.add_argument("gtf_file")
    args = parser.parse_args()    
    categorize_transcripts(args.ref_gtf_file, args.gtf_file)

if __name__ == '__main__':
    main()
