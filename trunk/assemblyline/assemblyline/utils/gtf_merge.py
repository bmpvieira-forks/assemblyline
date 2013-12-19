'''
Created on Dec 18, 2013

@author: mkiyer
'''
'''
Created on Dec 17, 2013

@author: mkiyer
'''
import os
import sys
import argparse
import logging
import collections
import itertools
import operator

from assemblyline.lib.gtf import GTFFeature
from assemblyline.lib.bx.intersection import Interval, IntervalTree
from assemblyline.lib.bx.cluster import ClusterTree

class Feature(object):    
    __slots__ = ('chrom', 'start', 'end', 'strand', 'exons', 'attrs')

    def __init__(self):
        self.chrom = None
        self.start = None
        self.end = None
        self.strand = None
        self.attrs = {}
        self.exons = []

    def to_gtf_features(self, source=None):
        if source is None:
            source = 'source'
        # transcript feature
        f = GTFFeature()
        f.seqid = self.chrom
        f.source = source
        f.feature_type = 'transcript'
        f.start = self.start
        f.end = self.end
        f.score = 1000.0
        f.strand = self.strand
        f.phase = '.'
        f.attrs = self.attrs.copy()
        features = [f]
        # exon features
        for i,e in enumerate(self.exons):
            start,end = e
            f = GTFFeature()
            f.seqid = self.chrom
            f.source = source
            f.feature_type = 'exon'
            f.start = start
            f.end = end
            f.score = 1000.0
            f.strand = self.strand
            f.phase = '.'
            f.attrs = self.attrs.copy()
            f.attrs["exon_number"] = i
            features.append(f)
        return features

    @staticmethod
    def from_genepred(line, attr_fields=None):        
        if line is None:
            return None
        line = line.strip()
        if line.startswith('#'):
            return None
        if line.startswith('track'):
            return None
        fields = line.split('\t')
        assert len(fields) >= 10
        #name,chrom,strand,txStart,txEnd,cdsStart,cdsEnd,exonCount,exonStarts,exonEnds,name
        f = Feature()
        f.attrs['transcript_id'] = fields[0]
        f.attrs['gene_id'] = fields[0]
        f.chrom = fields[1]
        f.strand = fields[2]
        f.start = int(fields[3])
        f.end = int(fields[4])
        # coding status --> transcript_type attribute
        cds_start = int(fields[5])
        cds_end = int(fields[6])
        coding = (cds_start != cds_end)
        if coding:
            f.attrs['transcript_type'] = 'protein_coding'
        else:
            f.attrs['transcript_type'] = 'non_coding'
        # exons
        exon_starts = map(int, fields[8].split(',')[:-1])
        exon_ends = map(int, fields[9].split(',')[:-1])
        f.exons = []
        for start,end in itertools.izip(exon_starts, exon_ends):
            f.exons.append((start,end))
        if f.strand == '-':
            f.exons.reverse()
        if attr_fields is not None:
            for i,a in enumerate(attr_fields):
                f.attrs[a] = fields[10+i]
        return f

    @staticmethod
    def parse_genepred(line_iter, attr_fields=None):
        rename_dict = collections.defaultdict(lambda: 0)
        for line in line_iter:
            if not line:
                continue
            if not line.strip():
                continue            
            if line.startswith("#"):
                continue
            if line.startswith("track"):
                continue
            f = Feature.from_genepred(line, attr_fields)
            t_id = f.attrs['transcript_id']
            if t_id in rename_dict:
                count = rename_dict[t_id]
                new_t_id = '%s.%d' % (t_id,count)
                f.attrs['transcript_id'] = new_t_id
                logging.warning('transcript id "%s" already present so renamed to "%s"' % (t_id, new_t_id))
            rename_dict[t_id] += 1
            yield f

    @staticmethod
    def from_gtf_features(features):
        self = Feature()
        self.attrs['transcript_id'] = features[0].attrs['transcript_id']
        self.attrs['gene_id'] = features[0].attrs['gene_id']
        self.attrs['gene_name'] = features[0].attrs['gene_name']
        self.attrs['transcript_type'] = features[0].attrs['transcript_type']
        self.attrs['source'] = features[0].source
        self.chrom = features[0].seqid
        self.strand = features[0].strand
        self.start = min(f.start for f in features)
        self.end = max(f.end for f in features)
        self.exons = []
        for f in features:
            self.exons.append((f.start, f.end))
        self.exons.sort()
        if self.strand == '-':
            self.exons.reverse()
        return self
        
    @staticmethod
    def parse_gtf(filename):
        # read all transcripts
        t_dict = collections.defaultdict(lambda: [])
        i = 0
        for f in GTFFeature.parse(open(filename)):
            i += 1
            if (i % 100000) == 0:
                logging.debug('parse_gtf read %d lines' % (i))
            if f.feature_type != 'exon':
                continue
            t_id = f.attrs['transcript_id']
            t_dict[t_id].append(f)
        i = 0
        for features in t_dict.itervalues():
            yield Feature.from_gtf_features(features)
            i += 1
        logging.debug('Parsed %d transcripts' % (i))

def cluster_gtf_features(gtf_files, source=None):
    # read all features
    chrom_feature_dict = collections.defaultdict(lambda: collections.defaultdict(lambda: []))    
    for gtf_file in gtf_files:
        logging.debug('Parsing gtf file: %s' % (gtf_file))
        for f in Feature.parse_gtf(gtf_file):
            # bin by chromosome and strand
            chrom_feature_dict[f.chrom][f.strand].append(f)
    # cluster transcripts into genes
    if source is None:
        source = 'merge'
    logging.debug('Clustering transcripts into genes')
    cur_gene_id = 1
    for strand_feature_dict in chrom_feature_dict.itervalues():
        for strand_features in strand_feature_dict.itervalues():
            # initialize each transcript to be in a 'gene' by itself
            cluster_map = {}
            cluster_tree = ClusterTree(0,1)
            for i,f in enumerate(strand_features):
                cluster_map[i] = set((i,))
                for start,end in f.exons:
                    cluster_tree.insert(start, end, i)
            for start, end, indexes in cluster_tree.getregions():
                # group transcripts into larger clusters
                new_cluster = set()
                for i in indexes:
                    new_cluster.update(cluster_map[i])
                # reassign transcript clusters to new cluster
                for i in indexes:
                    cluster_map[i] = new_cluster
            del cluster_tree
            # now all transcripts are assigned to a gene cluster
            # enumerate all gene clusters
            clusters = set()
            for clust in cluster_map.itervalues():
                clusters.add(frozenset(clust))
            del cluster_map
            # now assign gene ids to each cluster
            for clust in clusters:
                new_gene_id = '%s%011d' % (source, cur_gene_id)
                for i in clust:
                    f = strand_features[i]
                    f.attrs['orig_gene_id'] = f.attrs['gene_id']
                    f.attrs['gene_id'] = new_gene_id
                cur_gene_id += 1
    # output genes
    logging.debug('Writing transcripts')
    for chrom in sorted(chrom_feature_dict):
        strand_feature_dict = chrom_feature_dict[chrom]
        features = []
        for strand_features in strand_feature_dict.itervalues():
            features.extend(strand_features)
        features.sort(key=operator.attrgetter('start'))
        for f in features:
            for gtf_feature in f.to_gtf_features(source=source):
                yield str(gtf_feature)

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument('--source', dest='source', default=None)
    parser.add_argument('gtf_file_list')
    args = parser.parse_args()
    gtf_file_list = args.gtf_file_list
    gtf_files = []
    with open(gtf_file_list) as fileh:
        for line in fileh:
            filename = line.strip()
            if not filename:
                continue
            if not os.path.exists(filename):
                parser.error('GTF file "%s" not found' % (filename))
            gtf_files.append(filename)
            logging.info('gtf file: %s' % (filename))
    for f in cluster_gtf_features(gtf_files, source=args.source):
        print f
    logging.info("Done")
    return 0

if __name__ == '__main__':
    sys.exit(main())
