'''
Created on Feb 21, 2012

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
import sys
import os
import logging
import argparse
import numpy as np
import collections

from assemblyline.lib.bx.cluster import ClusterTree
from assemblyline.lib.base import Library
from assemblyline.lib.gtf import GTFFeature

htseq_special_fields = set(["no_feature", "ambiguous", "too_low_aQual", "not_aligned", "alignment_not_unique"])

def get_gtf_metadata(gtf_file, 
                      omit_attrs=None,
                      group_by="gene_id", 
                      feature_type="exon"):
    if omit_attrs is None:
        omit_attrs = []
    # read gtf file and group by gene
    gene_feature_map = collections.defaultdict(lambda: [])
    gene_attrs_set = set()
    for feature in GTFFeature.parse(open(gtf_file)):
        if feature.feature_type != feature_type:
            continue
        feature_id = feature.attrs[group_by]
        gene_feature_map[feature_id].append(feature)
        gene_attrs_set.update(feature.attrs.keys())
    gene_attrs_set.difference_update(omit_attrs)
    gene_attrs_list = sorted(gene_attrs_set)
    metadata_fields = ["tracking_id", "locus", "strand", "num_exons", "transcript_length"] + gene_attrs_list
    metadata_inds = dict((x,i) for i,x in enumerate(metadata_fields))
    metadata_dict = {}
    # output metadata sorted by gene id
    for feature_id,features in gene_feature_map.iteritems():
        # collect attributes for this gene
        attrdict = collections.defaultdict(lambda: set())
        # cluster exons together for each gene
        cluster_tree = ClusterTree(0,1)
        for i,f in enumerate(features):
            cluster_tree.insert(f.start, f.end, i)
            for k,v in f.attrs.iteritems():
                if k in gene_attrs_set:
                    # some attributes have multiple values separated by a comma
                    attrdict[k].update(v.split(','))
        # determine larger exon clusters
        transcript_length = 0
        exon_clusters = []
        for start, end, indexes in cluster_tree.getregions():
            exon_clusters.append((start,end))
            transcript_length += (end - start)
        del cluster_tree
        chrom = features[0].seqid
        locus_start = min(e[0] for e in exon_clusters)
        locus_end = max(e[1] for e in exon_clusters)
        locus_string = "%s:%d-%d" % (chrom, locus_start, locus_end)
        strand = features[0].strand
        num_exons = len(exon_clusters)
        # make metadata row
        metadata = [feature_id, locus_string, strand, num_exons, transcript_length] + ['NA'] * len(gene_attrs_list)
        # get all attributes
        for k,vals in attrdict.iteritems():
            ind = metadata_inds[k]
            metadata[ind] = ','.join(map(str, sorted(vals)))
        metadata_dict[metadata[0]] = metadata
    return metadata_fields, metadata_dict

def read_htseq_count_file(filename, genes, specials):
    gcounts = []
    scounts = []
    gene_ind = 0
    special_ind = 0
    for line in open(filename):
        fields = line.strip().split('\t')
        if fields[0] in htseq_special_fields:
            assert fields[0] == specials[special_ind]
            scounts.append(int(fields[1]))
            special_ind += 1
        else:
            assert fields[0] == genes[gene_ind]
            gcounts.append(int(fields[1]))
            gene_ind += 1
    return gcounts, scounts

def read_htseq_id_order(filename):
    genes = []
    specials = []
    for line in open(filename):
        fields = line.strip().split('\t')        
        if fields[0] in htseq_special_fields:
            specials.append(fields[0])
        else:
            genes.append(fields[0])
    return genes, specials

def main():
    # Command line parsing
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument('--omit', dest="omit_list",
                        default="exon_number,score,frac")
    parser.add_argument('--htseq-filename', dest="htseq_filename", default="htseq.txt")
    parser.add_argument('gtf_file')
    parser.add_argument('library_table')
    parser.add_argument('input_dir')
    parser.add_argument('output_file_prefix')
    args = parser.parse_args()
    if not os.path.exists(args.input_dir):
        parser.error("input dir %s not found" % (args.input_dir))
    if not os.path.exists(args.gtf_file):
        parser.error("gtf file %s not found" % (args.gtf_file))
    input_dir = args.input_dir
    htseq_filename = args.htseq_filename
    omit_list = args.omit_list.split(",")
    pheno_file = args.output_file_prefix + "_pheno_data.txt"
    count_file = args.output_file_prefix + "_count_data.txt"
    # aggregate runs into a single gene expression matrix
    logging.info("Reading library info and count data")
    count_vectors = []
    genes = None
    specialgenes = None
    library_ids = []
    infileh = open(args.library_table)    
    # first line is header fields
    header_fields = infileh.next().strip().split('\t')
    # figure out which column has library id
    lib_id_col = header_fields.index("library_id")
    phenofileh = open(pheno_file, "w")
    for line in infileh:
        fields = line.strip().split('\t')
        library_id = fields[lib_id_col]
        htseq_file = os.path.join(input_dir, library_id, htseq_filename)
        if not os.path.exists(htseq_file):
            logging.warning("[SKIPPED] Library %s missing count file" % (library_id))
            continue
        # get the order of the genes in the first file
        if genes is None:
            genes, specialgenes = read_htseq_id_order(htseq_file)
            # output pheno file header
            pheno_header_fields = list(header_fields)
            pheno_header_fields.extend(specialgenes)
            pheno_header_fields.append("total_counts")
            print >>phenofileh, '\t'.join(pheno_header_fields)
        # read counts
        logging.info("Processing library %s" % (library_id))
        counts, specialcounts = read_htseq_count_file(htseq_file, genes, specialgenes)
        count_vectors.append(counts)
        library_ids.append(library_id)
        # write pheno file line
        pheno_fields = list(fields)
        pheno_fields.extend(specialcounts)
        pheno_fields.append(sum(counts))
        print >>phenofileh, '\t'.join(map(str, pheno_fields))
    infileh.close()
    phenofileh.close()
    if len(count_vectors) == 0:
        logging.error("No valid libraries, exiting")
        return 1
    # combine frag count vectors to get expression matrix
    logging.info("Merging count vectors")
    mat = np.vstack(count_vectors).T
    del count_vectors
    # get gene metadata
    logging.info("Reading gene metadata")
    metadata_fields, metadata_dict = get_gtf_metadata(args.gtf_file, omit_attrs=omit_list)
    # write matrix file
    logging.info("Writing count file")
    header_fields = list(metadata_fields)
    header_fields.extend(library_ids)
    fileh = open(count_file, "w")
    print >>fileh, '\t'.join(header_fields)
    for i,feature_id in enumerate(genes):
        fields = list(map(str, metadata_dict[feature_id]))
        fields.extend(map(str, mat[i,:]))
        print >>fileh, '\t'.join(fields)
    fileh.close()
    return 0

if __name__ == '__main__':
    sys.exit(main())
