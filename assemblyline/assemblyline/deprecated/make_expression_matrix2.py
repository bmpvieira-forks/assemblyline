'''
Created on Sep 17, 2013

@author: mkiyer
'''
import sys
import os
import logging
import argparse
import collections
import re
import numpy as np

from assemblyline.lib.gtf import GTFFeature

#COUNTS = FPKM*transcript_length/1000000000 * MASS

MAP_MASS_RE = re.compile(r'Normalized Map Mass: (.+)$')
FPKMTuple = collections.namedtuple('FPKMTuple', 
                                   ('coverage', 
                                    'FPKM', 
                                    'FPKM_conf_lo', 
                                    'FPKM_conf_hi', 
                                    'FPKM_status'))
GTF_ATTRS = ['transcript_id',
             'gene_id',
             'tss_id',
             'category', 
             'nearest_dist',
             'annotation_sources',
             'nearest_gene_names',
             'nearest_gene_ids']
GTF_ATTRS_SET = set(GTF_ATTRS)

class TranscriptMetadata(object):
    def __init__(self):
        self.chrom = None
        self.start = 0
        self.end = 0
        self.strand = '.'
        self.num_exons = 0
        self.length = 0
        for attr in GTF_ATTRS:
            setattr(self, attr, '')

def get_gtf_metadata(gtf_file, metadata_file):
    # read gtf file and group by gene
    metadata_dict = {}
    for feature in GTFFeature.parse(open(gtf_file)):
        if feature.feature_type != "exon":
            continue
        t_id = feature.attrs["transcript_id"]
        if t_id not in metadata_dict:
            # instantiate new metadata
            m = TranscriptMetadata()
            m.chrom = feature.seqid
            m.strand = feature.strand
            m.start = feature.start
            m.end = feature.end
            for attr in GTF_ATTRS_SET:
                setattr(m, attr, feature.attrs.get(attr, ''))
            metadata_dict[t_id] = m
        else:
            m = metadata_dict[t_id]
        # update metadata
        m.start = feature.start if feature.start < m.start else m.start
        m.end = feature.end if feature.end > m.end else m.end
        m.length += (feature.end - feature.start)
        m.num_exons += 1
    fileh = open(metadata_file, 'w')
    header_fields = ['tracking_id', 'locus', 'strand', 'num_exons', 'transcript_length'] + GTF_ATTRS
    print >>fileh, '\t'.join(header_fields)
    tracking_ids = sorted(metadata_dict)
    for t_id in tracking_ids:
        m = metadata_dict[t_id]
        fields = [t_id,
                  '%s:%d-%d' % (m.chrom, m.start, m.end),
                  m.strand,
                  m.num_exons,
                  m.length]
        for attr in GTF_ATTRS:
            fields.append(getattr(m, attr))
        print >>fileh, '\t'.join(map(str, fields))
    fileh.close()
    return tracking_ids

def get_cufflinks_map_mass(cufflinks_log_file):
    for line in open(cufflinks_log_file):
        if not line.startswith('>'):
            continue
        m = MAP_MASS_RE.search(line)
        if not m:
            continue
        map_mass = float(m.group(1))
        return map_mass

def get_library_metadata(library_table_file, input_dir, pheno_file):
    infileh = open(library_table_file)    
    # first line is header fields
    header_fields = infileh.next().strip().split('\t')
    # open pheno file
    phenofileh = open(pheno_file, "w")
    header_fields.append("cufflinks_map_mass")
    print >>phenofileh, '\t'.join(header_fields)
    # figure out which column has library id
    library_id_col = header_fields.index("library_id")
    library_ids = []
    for line in infileh:
        if line.startswith("#"):
            continue
        if not line:
            continue
        line = line.strip()
        if not line:
            continue
        fields = line.strip().split('\t')
        library_id = fields[library_id_col]
        # get normalized map mass
        cufflinks_done_file = os.path.join(input_dir, library_id, "cufflinks.done")
        if not os.path.exists(cufflinks_done_file):
            logging.warning("[SKIPPED] Library %s missing 'cufflinks.done' file" % (library_id))
            continue
        cufflinks_log_file = os.path.join(input_dir, library_id, "cufflinks.log")
        map_mass = get_cufflinks_map_mass(cufflinks_log_file)
        library_ids.append(library_id)
        # write pheno file line
        pheno_fields = list(fields)
        pheno_fields.append(str(map_mass))
        print >>phenofileh, '\t'.join(pheno_fields)
    phenofileh.close()
    infileh.close()
    return library_ids

def get_isoform_fpkm_data(filename):
    fileh = open(filename)
    header_fields = fileh.next().strip().split('\t')
    tracking_id_ind = header_fields.index('tracking_id')
    fpkm_ind = header_fields.index('FPKM')
    fpkm_status_ind = header_fields.index('FPKM_status')
    data_dict = {}
    for line in fileh:
        fields = line.strip().split('\t')
        tracking_id = fields[tracking_id_ind]        
        fpkm = float(fields[fpkm_ind])
        status = fields[fpkm_status_ind]
        if status != "OK":
            fpkm = -1.0
        data_dict[tracking_id] = fpkm
    fileh.close()
    return data_dict

def main():
    # Command line parsing
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", dest="mode", choices=['htseq', 'cufflinks'],
                        default='htseq')
    parser.add_argument('--omit', dest="omit_list",
                        default="exon_number,score,frac")
    parser.add_argument('gtf_file')
    parser.add_argument('library_table')
    parser.add_argument('input_dir')
    parser.add_argument('output_dir')
    args = parser.parse_args()
    if not os.path.exists(args.gtf_file):
        parser.error("GTF file '%s' not found" % (args.gtf_file))
    if not os.path.exists(args.library_table):
        parser.error("Library table file '%s' not found" % (args.library_table))
    if not os.path.exists(args.input_dir):
        parser.error("Input directory '%s' not found" % (args.input_dir))
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
    pheno_file = os.path.join(args.output_dir, 'phenos.txt')
    metadata_file = os.path.join(args.output_dir, 'metadata.txt')
    # collect transcript metadata from GTF file
    logging.info("Reading GTF file")
    tracking_ids = get_gtf_metadata(args.gtf_file, metadata_file)
    logging.info("\tfound %d transcripts" % (len(tracking_ids)))
    # collect library ids and write phenotypes file
    logging.info("Reading library table file")
    library_ids = get_library_metadata(args.library_table, args.input_dir, pheno_file)
    logging.info("\tfound %d libraries" % (len(library_ids)))
    # create matrix
    logging.info("Creating expression matrix")
    mat = np.zeros((len(tracking_ids),len(library_ids)), dtype='float32')
    for j,library_id in enumerate(library_ids):
        logging.debug("\tlibrary %d/%d: '%s'" % (j+1,len(library_ids),library_id))
        isoform_fpkm_file = os.path.join(args.input_dir, library_id, 'isoforms.fpkm_tracking')
        fpkm = get_isoform_fpkm_data(isoform_fpkm_file)
        fpkm = np.array([fpkm.get(x,-2.0) for x in tracking_ids], dtype='float32')
        mat[:,j] = fpkm
    # write to memmap
    matrix_file = os.path.join(args.output_dir, 'isoform_fpkm.mmap')
    fp = np.memmap(matrix_file, dtype='float32', mode='w+', 
                   shape=(len(tracking_ids),len(library_ids)))
    fp[:] = mat[:]
    del fp
    # write to file
    fileh = open(os.path.join(args.output_dir, 'isoform_fpkm.txt'), 'w')
    header_fields = ['tracking_id']
    header_fields.extend(library_ids)
    print >>fileh, '\t'.join(header_fields)    
    def fpkm2str(x):
        return "NA" if x < 0 else str(x)
    for i,tracking_id in enumerate(tracking_ids):
        fields = [tracking_id]
        fields.extend(map(fpkm2str, mat[i,:]))
        print >>fileh, '\t'.join(fields)
    fileh.close()
    logging.info("Done.")
    return 0


if __name__ == '__main__':
    sys.exit(main())
