'''
Created on Dec 11, 2011

@author: mkiyer
'''
import collections
import sys
import os
import logging
import argparse
from multiprocessing import Pool
import numpy as np

from assemblyline.lib.bx.cluster import ClusterTree
from assemblyline.lib.sampletable import LibInfo
from assemblyline.lib.quantify import quantify_gene_abundance
from assemblyline.lib.transcript_parser import parse_gtf
from assemblyline.lib.transcript import Exon, POS_STRAND, NEG_STRAND, NO_STRAND, strand_str_to_int, strand_int_to_str
from assemblyline.lib import genome

def run_quantify(bam_file, gtf_file, output_dir, sample_id, library_id):
    output_file = os.path.join(output_dir, "gene_abundance.txt")
    # ensure BAM index can be found by pysam module    
    bai_file = bam_file + ".bai"
    split_bai_file = os.path.splitext(bam_file)[0] + ".bai"
    made_symlink = False
    if not os.path.exists(bai_file):
        if not os.path.exists(split_bai_file):
            logging.error("BAM index not found")
            return 1, sample_id, library_id
        os.symlink(split_bai_file, bai_file)
        made_symlink = True
    # run abundance estimation
    try:
        quantify_gene_abundance(gtf_file, bam_file, output_file)
    except Exception as e:
        logging.error("Exception caught: %s" % (str(e)))
        return 2, sample_id, library_id
    if made_symlink:
        # os.remove(bai_file)
        pass
    # create job.done file
    open(os.path.join(output_dir, "job.done"), "w").close()
    return 0, sample_id, library_id

def run_quantify_task(args):
    return run_quantify(**args)

#def read_gene_metadata(filename):
#    f = open(filename)
#    # skip header
#    f.next()
#    gene_metadata  = []    
#    for line in f:
#        if line.startswith("#"):
#            continue
#        fields = line.strip().split('\t')
#        gene_id = fields[0]
#        strand = fields[1]
#        chrom = fields[2]
#        chrom = genome.ensembl_to_ucsc(chrom)
#        locus = "%s:%s-%s" % (chrom, fields[3], fields[4])
#        nearest_ref_id = "-"
#        class_code = "strand=%s,num_exons=%s" % (strand, fields[5])
#        transcript_length = fields[6]
#        gene_metadata.append((gene_id, locus, nearest_ref_id, class_code, transcript_length))
#    f.close()
#    return gene_metadata



def get_gene_metadata(locus_chrom, locus_start, locus_end, transcripts):
    class Gene(object):
        def __init__(self):
            self.strand = NO_STRAND
            self.exons = set()
    gene_id_map = collections.defaultdict(lambda: Gene())
    for t in transcripts:
        g = gene_id_map[t.gene_id]
        g.exons.update(t.exons)
        g.strand = t.strand
    # cluster exons together for each gene
    for g in gene_id_map.itervalues():
        cluster_tree = ClusterTree(0,1)
        for i,e in enumerate(g.exons):
            cluster_tree.insert(e.start, e.end, i)
        exon_clusters = []
        for start, end, indexes in cluster_tree.getregions():
            exon_clusters.append(Exon(start,end))
        # update exons
        g.exons = exon_clusters
        del cluster_tree
    # get transcript attributes
    category = transcripts[0].attrs.get("category", "-")
    if not category:
        category = "-"
    nearest_genes = transcripts[0].attrs.get("nearest_genes", "-")
    if not nearest_genes:
        nearest_genes = "-"
    #annotation_sources = transcripts[0].attrs.get("annotation_sources", "")
    for gene_id, g in gene_id_map.iteritems():
        # sum mass on all strands over the exonic intervals of the transcript
        length = sum((e.end - e.start) for e in g.exons)
        strand = strand_int_to_str(g.strand)
        description = "category=%s,strand=%s,exons=%d" % (category, strand, len(g.exons))
        locus_string = "%s:%d-%d" % (locus_chrom, g.exons[0].start, g.exons[-1].end) 
        yield (gene_id, locus_string, nearest_genes, description, length)

def get_metadata(gtf_file):
    for transcripts in parse_gtf(open(gtf_file)):
        chrom = transcripts[0].chrom
        start = transcripts[0].start
        end = max(tx.end for tx in transcripts)
        for metadata in get_gene_metadata(chrom, start, end, transcripts):
            yield metadata

def read_gene_abundance_vector(filename, gene_metadata):
    f = open(filename)
    # skip header
    f.next()
    frags = []
    line_num = 0
    for line in f:
        if line.startswith("#"):
            continue
        fields = line.strip().split('\t')
        assert fields[0] == gene_metadata[line_num][0]
        line_num += 1
        frags.append(int(fields[7]))
    f.close()
    return np.array(frags, dtype=int)

def main():
    # Command line parsing
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--output-dir', dest="output_dir", default=".")
    parser.add_argument('-p', '--num-processors', dest="num_processors", type=int, default=1)
    parser.add_argument('gtf_file')
    parser.add_argument('library_table')
    parser.add_argument('matrix_file')
    args = parser.parse_args()
    # check args
    num_processes = max(1, args.num_processors)
    gtf_file = os.path.abspath(args.gtf_file)
    # read gene information
    logging.info("Reading gene metadata")
    gene_metadata = list(get_metadata(gtf_file))
    # parse library table
    logging.info("Parsing library table")
    tasks = []
    libinfos = []
    for libinfo in LibInfo.from_file(args.library_table):
        if not libinfo.is_valid():
            logging.warning("\tskipping cohort=%s patient=%s sample=%s lib=%s lanes=%s" % 
                            (libinfo.cohort, libinfo.patient, libinfo.sample, 
                             libinfo.library, libinfo.lanes))
            continue
        libinfos.append(libinfo)
        lib_dir = os.path.join(args.output_dir, libinfo.sample, libinfo.library)
        if not os.path.exists(lib_dir):
            logging.debug("Creating directory %s" % (lib_dir))
            os.makedirs(lib_dir)
        # check for output
        job_done_file = os.path.join(lib_dir, "job.done")
        msg = "Adding task for sample=%s library=%s" % (libinfo.sample, libinfo.library)
        if os.path.exists(job_done_file):
            logging.info("[SKIPPED] %s" % msg)
        else:
            logging.info(msg)
            arg_dict = {'bam_file': libinfo.bam_file,
                        'gtf_file': gtf_file,
                        'output_dir': lib_dir,
                        'library_id': libinfo.library,
                        'sample_id': libinfo.sample}
            tasks.append(arg_dict)
    # use multiprocessing to parallelize
    logging.info("Running tasks")
    pool = Pool(processes=num_processes)
    result_iter = pool.imap_unordered(run_quantify_task, tasks)
    for retcode,sample_id,library_id in result_iter:        
        logging.debug("\tfinished sample %s library %s with return code %d" % (sample_id, library_id, retcode))
    pool.close()
    # aggregate runs into a single gene expression matrix
    # (assume this will fit into memory)
    logging.info("Building expression data matrix")
    frag_vectors = []
    for libinfo in libinfos:
        print libinfo.sample, libinfo.library
        lib_dir = os.path.join(args.output_dir, libinfo.sample, libinfo.library)
        gene_abundance_file = os.path.join(lib_dir, "gene_abundance.txt")
        job_done_file = os.path.join(lib_dir, "job.done")
        if not os.path.exists(job_done_file):
            logging.warning("skipping %s %s" % (libinfo.sample, libinfo.library))
            continue
        frag_vectors.append(read_gene_abundance_vector(gene_abundance_file, gene_metadata))    
    # combine frag count vectors to get expression matrix
    mat = np.vstack(frag_vectors).T
    del frag_vectors
    # write matrix file
    logging.info("Writing output file")
    f = open(args.matrix_file, "w")
    header_fields = ["tracking_id", "locus", "nearest_ref_id", "class_code", "transcript_length"]
    num_metadata_rows = len(header_fields)
    header_fields.extend([libinfo.library for libinfo in libinfos])
    print >>f, '\t'.join(header_fields)
    # sample names
    fields = [""] * (num_metadata_rows - 1)
    fields.append("sample")
    fields.extend([libinfo.sample for libinfo in libinfos])
    print >>f, '\t'.join(fields)    
    # sum of fragment vectors
    fields = [""] * (num_metadata_rows - 1)
    fields.append("total_frags")
    fields.extend(map(str,np.sum(mat,axis=0)))
    print >>f, '\t'.join(fields)    
    for i in xrange(len(gene_metadata)):
        fields = []
        fields.extend(map(str, gene_metadata[i]))
        fields.extend(map(str, mat[i,:]))
        print >>f, '\t'.join(fields)
    f.close()
    

if __name__ == '__main__':
    sys.exit(main())