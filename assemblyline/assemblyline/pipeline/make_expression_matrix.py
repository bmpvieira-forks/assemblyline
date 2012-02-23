'''
Created on Feb 21, 2012

@author: mkiyer
'''
import sys
import os
import logging
import argparse
import numpy as np
import collections

from assemblyline.lib.bx.cluster import ClusterTree
from assemblyline.lib.sampletable import LibInfo
from assemblyline.lib.gtf import GTFFeature

chrom_names = [["chr1", "1"],
               ["chr2", "2"],
               ["chr3", "3"],
               ["chr4", "4"],
               ["chr5", "5"],
               ["chr6", "6"],
               ["chr7", "7"],
               ["chr8", "8"],
               ["chr9", "9"],
               ["chr10", "10"],
               ["chr11", "11"],
               ["chr12", "12"],
               ["chr13", "13"],
               ["chr14", "14"],
               ["chr15", "15"],
               ["chr16", "16"],
               ["chr17", "17"],
               ["chr18", "18"],
               ["chr19", "19"],
               ["chr20", "20"],
               ["chr21", "21"],
               ["chr22", "22"],
               ["chrX", "X"],
               ["chrY", "Y"],
               ["chrM", "MT"]]

ucsc_to_ensembl_dict = dict(chrom_names)
ensembl_to_ucsc_dict = dict((v,k) for k,v in ucsc_to_ensembl_dict.iteritems())

htseq_special_fields = set(["no_feature", "ambiguous", "too_low_aQual", "not_aligned", "alignment_not_unique"])

def get_gene_metadata(gtf_file, gene_id_attr="gene_id", feature_type="exon",
                      genome="ensembl"):
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
        if genome == "ensembl":
            chrom = ensembl_to_ucsc_dict.get(chrom, chrom)
        strand = features[0].strand
        start = min(x.start for x in features)
        end = max(x.end for x in features)
        gene_name = features[0].attrs['gene_name']
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
        locus_string = "%s:%d-%d" % (chrom, start, end)
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
    parser.add_argument('--cufflinks', dest="mode", action="store_const", const="cufflinks", default="htseq") 
    parser.add_argument('--htseq', dest="mode", action="store_const", const="htseq")
    parser.add_argument('--ucsc', dest="genome", action="store_const", const="ucsc", default="ensembl")
    parser.add_argument('--ensembl', dest="genome", action="store_const", const="ensembl")
    parser.add_argument('gtf_file')
    parser.add_argument('library_table')
    parser.add_argument('quantify_dir')
    parser.add_argument('output_file_prefix')
    args = parser.parse_args()
    # read library information
    logging.info("Reading library information")
    libinfos = []
    for libinfo in LibInfo.from_file(args.library_table):
        if not libinfo.is_valid():
            logging.warning("\tskipping cohort=%s patient=%s sample=%s lib=%s lanes=%s" % 
                            (libinfo.cohort, libinfo.patient, libinfo.sample, 
                             libinfo.library, libinfo.lanes))
            continue
        lib_dir = os.path.join(args.quantify_dir, libinfo.sample, libinfo.library)
        # check for output
        if args.mode == "htseq":
            count_file = os.path.join(lib_dir, "htseq_count.txt")
        elif args.mode == "cufflinks":
            count_file = os.path.join(lib_dir, "genes.expr")
        job_done_file = os.path.join(lib_dir, "job.done")
        if (os.path.exists(job_done_file) and os.path.exists(count_file)):
            libinfos.append(libinfo)
        else:
            logging.info("[SKIPPED] %s %s" % (libinfo.sample, libinfo.library))
    # get gene metadata
    logging.info("Reading gene features")
    gene_metadata_dict = {}
    for gene_metadata in get_gene_metadata(args.gtf_file, genome=args.genome):
        gene_metadata_dict[gene_metadata[0]] = gene_metadata
    # aggregate runs into a single gene expression matrix
    # (assume this will fit into memory)
    logging.info("Building expression data matrix")
    count_vectors = []
    special_vectors = []
    for libinfo in libinfos:
        logging.debug("sample: %s lib: %s" % (libinfo.sample, libinfo.library))
        lib_dir = os.path.join(args.quantify_dir, libinfo.sample, libinfo.library)
        count_file = os.path.join(lib_dir, "htseq_count.txt")
        gene_count_map, special_count_map = read_htseq_count_file(count_file)
        counts = []
        for gene_id in sorted(gene_metadata_dict):
            counts.append(gene_count_map[gene_id])
        specials = []
        for special_id in sorted(special_count_map):
            specials.append(special_count_map[special_id])
        count_vectors.append(counts)
        special_vectors.append(specials)
    # write phenotype file
    logging.info("Writing phenotype file")
    pheno_file = args.output_file_prefix + "_pheno_data.txt"
    f = open(pheno_file, "w")
    header_fields = ["cohort", "patient", "sample", "library", "lanes", "library_type", 
                     "has_pe_lanes", "frag_len_mean", "frag_len_std_dev"]
    header_fields.extend(sorted(htseq_special_fields))
    header_fields.append("total_counts")
    print >>f, '\t'.join(map(str, header_fields))
    for i,libinfo in enumerate(libinfos):
        fields = [libinfo.cohort, libinfo.patient, libinfo.sample, libinfo.library,
                  libinfo.lanes, libinfo.library_type, libinfo.has_pe_lanes, 
                  libinfo.frag_len_mean, libinfo.frag_len_std_dev]
        fields.extend(special_vectors[i])
        fields.append(sum(count_vectors[i]))
        print >>f, '\t'.join(map(str, fields))
    f.close()
    # write matrix file
    logging.info("Writing output file")
    # combine frag count vectors to get expression matrix
    mat = np.vstack(count_vectors).T
    del count_vectors
    matrix_file = args.output_file_prefix + "_count_data.txt"
    f = open(matrix_file, "w")
    header_fields = ["tracking_id", "locus", "nearest_ref_id", "class_code", "transcript_length"]
    header_fields.extend([libinfo.library for libinfo in libinfos])
    print >>f, '\t'.join(header_fields)
    for i,gene_id in enumerate(sorted(gene_metadata_dict)):
        fields = []
        fields.extend(map(str, gene_metadata_dict[gene_id]))
        fields.extend(map(str, mat[i,:]))
        print >>f, '\t'.join(fields)
    f.close()
    # attributes
    #num_metadata_rows = len(header_fields)
    #for attrname in ('cohort', 'patient', 'sample', 'library', 'lanes', 
    #                 'library_type', 'frag_len_mean', 'frag_len_std_dev'):
    #    fields = [""] * (num_metadata_rows - 1)
    #    fields.append(attrname)
    #    fields.extend([getattr(libinfo, attrname) for libinfo in libinfos])
    #    print >>f, '\t'.join(fields)
    # special attributes
    #if args.mode == "htseq":
    #    special_mat = np.vstack(special_vectors).T
    #    for i,attrname in enumerate(sorted(htseq_special_fields)):
    #        fields = [""] * (num_metadata_rows - 1)
    #        fields.append(attrname)
    #        fields.extend(map(str, special_mat[i,:]))
    #        print >>f, '\t'.join(fields)
    # sum of fragment vectors
    #fields = [""] * (num_metadata_rows - 1)
    #fields.append("total_frags")
    #fields.extend(map(str,np.sum(mat,axis=0)))
    #print >>f, '\t'.join(fields)
  

if __name__ == '__main__':
    sys.exit(main())
