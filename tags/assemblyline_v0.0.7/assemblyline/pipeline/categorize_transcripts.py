'''
Created on Jan 25, 2012

@author: mkiyer
'''
import logging
import argparse
import collections
import operator
import os

from assemblyline.lib.bx.intersection import Interval, IntervalTree
from assemblyline.lib.bx.cluster import ClusterTree

from assemblyline.lib import gtf
from assemblyline.lib.transcript import strand_int_to_str
from assemblyline.lib.transcript_parser import parse_gtf, cufflinks_attr_defs

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
    gene_map = collections.defaultdict(lambda: Gene())
    logging.info("Reading GTF file")
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
    return locus_trees

def categorize_gene_transcripts(transcripts, locus_trees):
    any_is_sense = False
    any_is_coding = False
    any_exon_hits = False
    any_transcript_hits = False
    gene_names = set()
    annotation_sources = set()
    for t in transcripts:
        # intersect transcript with reference loci
        transcript_hits = locus_trees[t.chrom].find(t.start, t.end)
        strand = strand_int_to_str(t.strand)
        if len(transcript_hits) > 0:
            # this transcript overlaps known loci, now categories
            # as sense/antisense, coding/noncoding, exon/intron
            for transcript_hit in transcript_hits:
                gene_tree = transcript_hit.value
                for exon in t.exons:
                    for gene_hit in gene_tree.find(exon.start, exon.end):
                        g = gene_hit.value
                        any_exon_hits = True
                        if cmp_strand(g.strand, strand):
                            any_is_sense = True
                        if g.is_coding:
                            any_is_coding = True
                        gene_names.update(g.gene_names)
                        annotation_sources.update(g.annotation_sources)
            any_transcript_hits = True
    if any_transcript_hits:        
        if any_exon_hits:
            if any_is_sense:
                if any_is_coding:
                    category = "protein"
                else:
                    category = "ncrna"
            else:
                category = "antisense"
        else:
            category = "intronic"
        nearest_genes = "|".join(sorted(gene_names))
        annotation_sources = "|".join(sorted(annotation_sources))
    else:
        category = "intergenic"
        gene_start = min(t.start for t in transcripts)
        gene_end = max(t.end for t in transcripts) 
        nearest_hit, nearest_genes, nearest_dist = None, "NA", 10000000
        before_hits = locus_trees[t.chrom].before(gene_start, num_intervals=1, max_dist=10000000)        
        for locus_hit in before_hits:
            dist = gene_start - locus_hit.end
            if dist < nearest_dist:
                nearest_dist = dist
                nearest_hit = locus_hit
        after_hits = locus_trees[t.chrom].after(gene_end, num_intervals=1, max_dist=10000000)
        for locus_hit in after_hits:
            dist = locus_hit.start - gene_end
            if dist < nearest_dist:
                nearest_dist = dist
                nearest_hit = locus_hit
        if nearest_hit is not None:        
            gene_tree = nearest_hit.value
            gene_names = set()
            for gene_hit in gene_tree.find(nearest_hit.start, nearest_hit.end):
                g = gene_hit.value
                gene_names.update(g.gene_names)
            nearest_genes = "|".join(sorted(gene_names)) 
        annotation_sources = "NA"
    for t in transcripts:
        t.attrs["category"] = category
        t.attrs["nearest_genes"] = nearest_genes
        t.attrs["annotation_sources"] = annotation_sources
        

def categorize_transcripts(ref_gtf_file, gtf_file):
    locus_trees = read_reference_gtf(ref_gtf_file)
    for locus_transcripts in parse_gtf(open(gtf_file), attr_defs=cufflinks_attr_defs):
        # group transcripts by gene id
        gene_transcript_map = collections.defaultdict(lambda: [])
        for t in locus_transcripts:
            gene_transcript_map[t.gene_id].append(t)
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
