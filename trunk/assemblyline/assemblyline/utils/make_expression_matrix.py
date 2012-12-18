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
from assemblyline.lib.librarytable import LibraryInfo
from assemblyline.lib.gtf import GTFFeature

htseq_special_fields = set(["no_feature", "ambiguous", "too_low_aQual", "not_aligned", "alignment_not_unique"])

def get_gene_metadata(gtf_file, 
                      gene_id_attr="gene_id", 
                      feature_type="exon"):
    # read gtf file and group by gene
    gene_feature_map = collections.defaultdict(lambda: [])
    for feature in GTFFeature.parse(open(gtf_file)):
        if feature.feature_type != feature_type:
            continue
        gene_id = feature.attrs[gene_id_attr]
        gene_feature_map[gene_id].append(feature)
    # output metadata sorted by gene id
    for gene_id in sorted(gene_feature_map):
        features = gene_feature_map[gene_id]
        chrom = features[0].seqid
        strand = features[0].strand
        start = min(x.start for x in features)
        end = max(x.end for x in features)
        if 'gene_name' in features[0].attrs:        
            gene_name = features[0].attrs['gene_name']
        else:
            gene_name = features[0].attrs['gene_id']
        # cluster exons together for each gene
        cluster_tree = ClusterTree(0,1)
        for i,f in enumerate(features):
            cluster_tree.insert(f.start, f.end, i)
        exon_clusters = []
        transcript_length = 0
        for start, end, indexes in cluster_tree.getregions():
            exon_clusters.append((start,end))
            transcript_length += (end - start)
        del cluster_tree
        num_exons = len(exon_clusters)
        locus_start = min(e[0] for e in exon_clusters)
        locus_end = max(e[1] for e in exon_clusters)
        locus_string = "%s:%d-%d" % (chrom, locus_start, locus_end)
        description = "category=%s,strand=%s,exons=%s" % (features[0].source, strand, num_exons)
        yield (gene_id, locus_string, gene_name, description, transcript_length)

def read_htseq_count_file(filename):
    gene_count_map = {}
    special_count_map = {}
    for line in open(filename):
        fields = line.strip().split('\t')
        if fields[0] in htseq_special_fields:
            special_count_map[fields[0]] = int(fields[1])
        else:
            gene_count_map[fields[0]] = int(fields[1])
    return gene_count_map, special_count_map

def main():
    # Command line parsing
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument('--input-filename', dest="input_filename", default="htseq.txt")
    parser.add_argument('gtf_file')
    parser.add_argument('library_table')
    parser.add_argument('input_dir')
    parser.add_argument('output_file_prefix')
    args = parser.parse_args()
    # get gene metadata
    logging.info("Reading gene features")
    gene_metadata_dict = {}
    for gene_metadata in get_gene_metadata(args.gtf_file):
        gene_metadata_dict[gene_metadata[0]] = gene_metadata
    # aggregate runs into a single gene expression matrix
    # (assume this will fit into memory)
    logging.info("Reading count data for each library")
    count_vector_dict = {}
    special_vector_dict = {}
    for lib in LibraryInfo.from_file(args.library_table):
        if not lib.is_valid():
            logging.info("[SKIPPED] Library %s invalid" % (lib.library_id))
            continue
        # check for input file
        lib_dir = os.path.join(args.input_dir, lib.library_id)
        count_file = os.path.join(lib_dir, args.input_filename)
        if not os.path.exists(count_file):
            logging.info("[SKIPPED] Library %s missing count file" % (lib.library_id))
            continue
        gene_count_map, special_count_map = read_htseq_count_file(count_file)
        counts = []
        for gene_id in sorted(gene_metadata_dict):
            counts.append(gene_count_map[gene_id])
        specials = []
        for special_id in sorted(special_count_map):
            specials.append(special_count_map[special_id])
        count_vector_dict[lib.library_id] = counts
        special_vector_dict[lib.library_id] = specials
    # write phenotype file
    logging.info("Writing phenotype file")
    fh = open(args.library_table)
    header_fields = fh.next().strip().split('\t')
    lib_id_col = header_fields.index("library_id")
    header_fields.extend(sorted(htseq_special_fields))
    header_fields.append("total_counts")
    pheno_file = args.output_file_prefix + "_pheno_data.txt"
    outfh = open(pheno_file, "w")
    print >>outfh, '\t'.join(map(str, header_fields))
    for line in fh:
        fields = line.strip().split('\t')
        lib_id = fields[lib_id_col]
        if lib_id not in count_vector_dict:
            continue
        fields.extend(special_vector_dict[lib_id])
        fields.append(sum(count_vector_dict[lib_id]))
        print >>outfh, '\t'.join(map(str, fields))
    outfh.close()
    fh.close()
    # write matrix file
    logging.info("Writing output file")
    # combine frag count vectors to get expression matrix
    sorted_lib_ids = sorted(count_vector_dict)
    mat = np.vstack(count_vector_dict[lib_id] for lib_id in sorted_lib_ids).T
    matrix_file = args.output_file_prefix + "_count_data.txt"
    f = open(matrix_file, "w")
    header_fields = ["tracking_id", "locus", "nearest_ref_id", "class_code", "transcript_length"]
    header_fields.extend(sorted_lib_ids)
    print >>f, '\t'.join(header_fields)
    for i,gene_id in enumerate(sorted(gene_metadata_dict)):
        fields = []
        fields.extend(map(str, gene_metadata_dict[gene_id]))
        fields.extend(map(str, mat[i,:]))
        print >>f, '\t'.join(fields)
    f.close()  

if __name__ == '__main__':
    sys.exit(main())
