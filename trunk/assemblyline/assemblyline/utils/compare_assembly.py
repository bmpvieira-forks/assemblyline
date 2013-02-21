'''
Created on Oct 8, 2012

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
import argparse
import logging
import os
import sys
import operator
import collections
import networkx as nx
import numpy as np

import assemblyline
from assemblyline.lib.base import GTFAttr
from assemblyline.lib.gtf import GTFFeature, sort_gtf
from assemblyline.lib.transcript import parse_gtf, Exon, NO_STRAND, POS_STRAND, NEG_STRAND
from assemblyline.lib.assemble.transcript_graph import find_exon_boundaries, split_exon, get_transcript_node_map
from assemblyline.lib.assemble.base import NODE_LENGTH, STRAND_SCORE, TRANSCRIPT_IDS 

TEST_IDS = 'test_ids'
REF_IDS = 'ref_ids'
TEST_SCORE = 'test_score'
REF_SCORE = 'ref_score'
NODE_ID_ATTRS = [TEST_IDS, REF_IDS]
STRAND_IS_REF = 'r'
STRAND_IS_ASSEMBLY = 'a'


def add_node_undirected(G, n, t, **kwargs):
    """
    add node to undirected graph
    """
    if n not in G: 
        attr_dict = {TRANSCRIPT_IDS: set(),
                     NODE_LENGTH: (n.end - n.start),
                     STRAND_SCORE: np.zeros(3,float)} 
        G.add_node(n, attr_dict=attr_dict)
    nd = G.node[n]
    nd[TRANSCRIPT_IDS].add(t.attrs[GTFAttr.TRANSCRIPT_ID])
    nd[STRAND_SCORE][t.strand] += t.score

def create_undirected_transcript_graph(transcripts, 
                                       add_node_func=add_node_undirected,
                                       **kwargs):
    '''
    add all transcripts to a single undirected graph
    '''
    # find the intron domains of the transcripts
    boundaries = find_exon_boundaries(transcripts)
    # initialize transcript graph as undirected at first
    G = nx.Graph()
    # add transcripts
    for t in transcripts:
        # split exons that cross boundaries and to get the
        # nodes in the transcript path
        nodes = []
        for exon in t.exons:
            for start,end in split_exon(exon, boundaries):
                nodes.append(Exon(start, end))
        # add nodes/edges to graph
        u = nodes[0]
        add_node_func(G, u, t, **kwargs)
        for v in nodes[1:]:
            add_node_func(G, v, t, **kwargs)
            G.add_edge(u, v)
            u = v
    return G

class MatchInfo(object):
    def __init__(self):
        self.has_pattern = 0
        self.num_introns = 0
        self.overlap = 0

class MatchStats(object):
    def __init__(self):
        self.gene_id = None
        self.t_id = None
        self.is_ref = None
        self.score = None
        self.length = None
        self.num_introns = None
        self.match_t_id = None
        self.match_pattern = None
        self.match_introns = None
        self.match_overlap = None

    @staticmethod
    def header_fields():
        return ["gene_id", "transcript_id", "is_ref", "score", "length",
                "num_introns", "match_transcript_id", "match_pattern",
                "match_introns", "match_overlap"]

    def __str__(self):
        fields = [self.gene_id, self.t_id, self.is_ref, 
                  self.score, self.length, self.num_introns, 
                  self.match_t_id, self.match_pattern, 
                  self.match_introns, self.match_overlap]
        return '\t'.join(map(str,fields))


class GlobalStats(object):
    fields = ('introns_both', 'introns_ref_only', 'introns_assembly_only',
              'patterns_both', 'patterns_ref_only', 'patterns_assembly_only', 
              'strand_cov_both', 'strand_cov_ref_only', 
              'strand_cov_assembly_only',
              'unstranded_cov_assembly_only',
              'cov_both', 'cov_ref_only', 'cov_assembly_only')
    def __init__(self):
        self.introns_both = 0
        self.introns_ref_only = 0
        self.introns_assembly_only = 0
        self.patterns_both = 0
        self.patterns_ref_only = 0
        self.patterns_assembly_only = 0
        self.strand_cov_both = 0
        self.strand_cov_ref_only = 0
        self.strand_cov_assembly_only = 0
        self.unstranded_cov_assembly_only = 0
        self.cov_both = 0
        self.cov_ref_only = 0
        self.cov_assembly_only = 0
    def __add__(self, x):
        z = GlobalStats()
        for f in GlobalStats.fields:
            val = getattr(self, f) + getattr(x, f)
            setattr(z, f, val)
        return z
    def report(self):
        # print stats report
        introns_total = self.introns_both + self.introns_ref_only + self.introns_assembly_only
        patterns_total = self.patterns_both + self.patterns_ref_only + self.patterns_assembly_only
        strand_cov_total = self.strand_cov_both + self.strand_cov_ref_only + self.strand_cov_assembly_only
        cov_total = self.cov_both + self.cov_ref_only + self.cov_assembly_only
        lines = ["introns_total=%d" % (introns_total),
                 "introns_both=%d" % (self.introns_both),
                 "introns_ref_only=%d" % (self.introns_ref_only),
                 "introns_assembly_only=%d" % (self.introns_assembly_only),
                 "introns_precision=%f" % (self.introns_both / float(self.introns_both + self.introns_assembly_only)),
                 "introns_recall=%f" % (self.introns_both / float(self.introns_both + self.introns_ref_only)),
                 "patterns_total=%d" % (patterns_total),
                 "patterns_both=%d" % (self.patterns_both),
                 "patterns_ref_only=%d" % (self.patterns_ref_only),
                 "patterns_assembly_only=%d" % (self.patterns_assembly_only),
                 "patterns_precision=%f" % (self.patterns_both / float(self.patterns_both + self.patterns_assembly_only)),
                 "patterns_recall=%f" % (self.patterns_both / float(self.patterns_both + self.patterns_ref_only)),
                 "strand_cov_total=%d" % (strand_cov_total),
                 "strand_cov_both=%d" % (self.strand_cov_both),
                 "strand_cov_ref_only=%d" % (self.strand_cov_ref_only),
                 "strand_cov_assembly_only=%d" % (self.strand_cov_assembly_only),
                 "strand_cov_precision=%f" % (self.strand_cov_both / float(self.strand_cov_both + self.strand_cov_assembly_only)),
                 "strand_cov_recall=%f" % (self.strand_cov_both / float(self.strand_cov_both + self.strand_cov_ref_only)),
                 "unstranded_cov_assembly_only=%d" % (self.unstranded_cov_assembly_only),
                 "cov_total=%d" % (cov_total),
                 "cov_both=%d" % (self.cov_both),
                 "cov_ref_only=%d" % (self.cov_ref_only),
                 "cov_assembly_only=%d" % (self.cov_assembly_only),
                 "cov_precision=%f" % (self.cov_both / float(self.cov_both + self.cov_assembly_only)),
                 "cov_recall=%f" % (self.cov_both / float(self.cov_both + self.cov_ref_only))]
        return '\n'.join(lines)

def add_node_match(G, n, t, **kwargs):
    if n not in G: 
        attr_dict = {TEST_IDS: set(),
                     REF_IDS: set(),
                     TEST_SCORE: 0.0,
                     REF_SCORE: 0.0}
        G.add_node(n, attr_dict=attr_dict)
    is_ref = bool(int(t.attrs[GTFAttr.REF]))
    t_id = t.attrs[GTFAttr.TRANSCRIPT_ID]
    nd = G.node[n]
    if is_ref:
        nd[REF_IDS].add(t_id)
        nd[REF_SCORE] += t.score
    else:
        nd[TEST_IDS].add(t_id)
        nd[TEST_SCORE] += t.score

def find_best_match(nodes, introns, G, other_node_attr, other_introns, 
                    other_patterns):
    match_dict = collections.defaultdict(lambda: MatchInfo())
    # search and count intron matches
    mypattern = []
    for start,end in introns:
        for match_t_id in other_introns[(start,end)]:
            match_dict[match_t_id].num_introns += 1
        mypattern.append((start,end))
    mypattern = tuple(mypattern)
    # find pattern matches
    if mypattern in other_patterns:
        for match_t_id in other_patterns[mypattern]:
            match_dict[match_t_id].has_pattern = 1
    # find and quantify base overlap
    match_node_map = collections.defaultdict(lambda: set())    
    for n in nodes:
        nd = G.node[n]
        for match_t_id in nd[other_node_attr]:
            match_node_map[match_t_id].add(n)
    for match_t_id,match_nodes in match_node_map.iteritems():
        overlap = sum((n.end - n.start) for n in match_nodes)
        match_dict[match_t_id].overlap = overlap
    if len(match_dict) == 0:
        return (None, 0, 0, 0)
    match_dict = dict(((m.has_pattern, m.num_introns, m.overlap), t_id) 
                      for t_id,m in match_dict.iteritems())
    sorted_match_keys = sorted(match_dict.keys(), reverse=True)
    best_match_tuple = sorted_match_keys[0]
    return (match_dict[best_match_tuple],) + best_match_tuple

def compare_locus(all_transcripts):
    # divide transcripts by strand
    strand_transcript_dicts = [{}, {}, {}]    
    for t in all_transcripts:
        t_id = t.attrs[GTFAttr.TRANSCRIPT_ID]        
        strand_transcript_dicts[t.strand][t_id] = t
    # process each strand separately
    for strand in xrange(0,3):
        strand_transcripts = strand_transcript_dicts[strand].values()
        # get introns and splicing patterns
        introns = [collections.defaultdict(lambda: set()),
                   collections.defaultdict(lambda: set())]
        patterns = [collections.defaultdict(lambda: set()),
                    collections.defaultdict(lambda: set())]
        transcripts = [[], []]
        for t in strand_transcripts:
            t_id = t.attrs[GTFAttr.TRANSCRIPT_ID]
            ref = int(t.attrs[GTFAttr.REF])
            transcripts[ref].append(t)
            pattern = []
            for start,end in t.iterintrons():
                introns[ref][(start,end)].add(t_id)
                pattern.append((start,end))
            pattern = tuple(pattern)
            patterns[ref][tuple(pattern)].add(t_id)
        # make transcript graph
        G = create_undirected_transcript_graph(strand_transcripts,
                                               add_node_func=add_node_match)
        # map transcripts to graph nodes
        test_node_map = get_transcript_node_map(G, node_attr=TEST_IDS)
        ref_node_map = get_transcript_node_map(G, node_attr=REF_IDS)
        node_maps = [test_node_map, ref_node_map]
        # compile statistics for each transcript
        for ref in xrange(0,2):
            not_ref = (ref + 1) % 2
            node_map = node_maps[ref]
            other_node_attr = NODE_ID_ATTRS[not_ref]
            other_introns = introns[not_ref]
            other_patterns = patterns[not_ref]
            for t in transcripts[ref]:
                g_id = t.attrs[GTFAttr.GENE_ID]                
                t_id = t.attrs[GTFAttr.TRANSCRIPT_ID]
                mynodes = node_map[t_id]                
                myintrons = list(t.iterintrons())
                match_tuple = find_best_match(mynodes, myintrons, G, 
                                              other_node_attr,
                                              other_introns,
                                              other_patterns)
                match_t_id, match_pattern, match_introns, match_overlap = match_tuple
                m = MatchStats()
                m.gene_id = g_id
                m.t_id = t_id
                m.is_ref = t.attrs[GTFAttr.REF]
                m.score = t.score
                m.length = t.length
                m.num_introns = len(myintrons)
                m.match_t_id = match_t_id
                m.match_pattern = match_pattern
                m.match_introns = match_introns
                m.match_overlap = match_overlap
                yield m

def add_node_cov_compare(G, n, t, **kwargs):
    if n not in G: 
        attr_dict = {STRAND_IS_REF: [False, False, False],
                     STRAND_IS_ASSEMBLY: [False, False, False]}
        G.add_node(n, attr_dict=attr_dict)
    nd = G.node[n]
    is_ref = bool(int(t.attrs[GTFAttr.REF]))
    if is_ref:
        nd[STRAND_IS_REF][t.strand] = True
    else:
        nd[STRAND_IS_ASSEMBLY][t.strand] = True

def gather_global_stats(transcripts):
    class CompareData(object):
        __slots__ = ('is_ref', 'is_assembly')
        def __init__(self):
            self.is_ref = False
            self.is_assembly = False
    # handle splicing related statistics
    strand_intron_dict = {POS_STRAND: collections.defaultdict(lambda: CompareData()),
                          NEG_STRAND: collections.defaultdict(lambda: CompareData()),
                          NO_STRAND: collections.defaultdict(lambda: CompareData())}
    strand_splicing_pattern_dict = {POS_STRAND: collections.defaultdict(lambda: CompareData()),
                                    NEG_STRAND: collections.defaultdict(lambda: CompareData()),
                                    NO_STRAND: collections.defaultdict(lambda: CompareData())}                                    
    # score and divide into ref/test transcript lists
    for t in transcripts:
        is_ref = bool(int(t.attrs[GTFAttr.REF]))
        intron_dict = strand_intron_dict[t.strand]
        pattern_dict = strand_splicing_pattern_dict[t.strand]
        pattern = []
        for start,end in t.iterintrons():
            k = (start,end)
            d = intron_dict[k]
            if is_ref:
                d.is_ref = True
            else:
                d.is_assembly = True
            pattern.append(k)
        if len(pattern) > 0:
            d = pattern_dict[tuple(pattern)]
            if is_ref:
                d.is_ref = True
            else:
                d.is_assembly = True
    # compile statistics
    stats_obj = GlobalStats()
    for strand, intron_dict in strand_intron_dict.iteritems():
        for d in intron_dict.itervalues():
            if d.is_ref and d.is_assembly:
                stats_obj.introns_both += 1
            elif d.is_ref:
                stats_obj.introns_ref_only += 1
            elif d.is_assembly:
                stats_obj.introns_assembly_only += 1
    for strand, splicing_pattern_dict in strand_splicing_pattern_dict.iteritems():
        for d in splicing_pattern_dict.itervalues():
            if d.is_ref and d.is_assembly:
                stats_obj.patterns_both += 1
            elif d.is_ref:
                stats_obj.patterns_ref_only += 1
            elif d.is_assembly:
                stats_obj.patterns_assembly_only += 1
    # coverage comparison
    # add all transcripts to a graph
    G = create_undirected_transcript_graph(transcripts,
                                           add_node_func=add_node_cov_compare)
    for n,nd in G.nodes_iter(data=True):
        length = (n.end - n.start)
        strand_is_ref = nd[STRAND_IS_REF]
        strand_is_assembly = nd[STRAND_IS_ASSEMBLY]
        for strand in xrange(0,2):
            is_ref = strand_is_ref[strand]
            is_assembly = strand_is_assembly[strand]
            if is_ref and is_assembly:
                stats_obj.strand_cov_both += length
            elif is_ref:
                stats_obj.strand_cov_ref_only += length
            elif is_assembly:
                stats_obj.strand_cov_assembly_only += length
        if strand_is_assembly[NO_STRAND]:
            stats_obj.unstranded_cov_assembly_only += length
        any_is_ref = any(strand_is_ref)
        any_is_assembly = any(strand_is_assembly)
        if any_is_ref and any_is_assembly:
            stats_obj.cov_both += length
        elif any_is_ref:
            stats_obj.cov_ref_only += length
        elif any_is_assembly:
            stats_obj.cov_assembly_only += length
    return stats_obj

def _parse_gtf_by_chrom(gtf_file):
    current_chrom = None
    exon_dict = collections.defaultdict(lambda: [])
    transcript_dict = {}
    for feature in GTFFeature.parse(open(gtf_file)):
        if (current_chrom != feature.seqid):
            if len(exon_dict) > 0:
                yield current_chrom, transcript_dict, exon_dict
                exon_dict = collections.defaultdict(lambda: [])
                transcript_dict = {}
            current_chrom = feature.seqid
        t_id = feature.attrs[GTFAttr.TRANSCRIPT_ID]
        if feature.feature_type == "transcript":
            transcript_dict[t_id] = feature
        elif feature.feature_type == "exon":
            exon_dict[t_id].append(feature)
    if len(exon_dict) > 0:
        yield current_chrom, transcript_dict, exon_dict

def add_gtf_file(gtf_file, outfh, is_ref, sample_id=None):
    refval = '1' if is_ref else '0'
    for chrom, transcript_dict, exon_dict in _parse_gtf_by_chrom(gtf_file):
        logging.debug("\tfinished chrom %s %d features" % (chrom, len(exon_dict)))
        # output reference transcripts
        for t_id, features in exon_dict.iteritems():
            # sort features (exons) by start position
            features.sort(key=operator.attrgetter('start'))
            # annotate exons as reference features
            for f in features:
                f.attrs[GTFAttr.REF] = refval
                print >>outfh, str(f)
            # transcript feature
            if t_id in transcript_dict:
                f = transcript_dict[t_id]
            else:
                f = GTFFeature()
                f.seqid = features[0].seqid
                f.source = features[0].source
                f.feature_type = 'transcript'
                f.start = features[0].start
                f.end = features[-1].end
                f.score = features[0].score
                f.strand = features[0].strand
                f.phase = '.'
                f.attrs = features[0].attrs.copy()
                if "exon_number" in f.attrs:
                    del f.attrs["exon_number"]
            f.attrs[GTFAttr.REF] = refval
            if sample_id is not None:
                f.attrs[GTFAttr.SAMPLE_ID] = sample_id
            print >>outfh, str(f)

def compare_assembly(ref_gtf_file, test_gtf_file, output_dir, 
                     gtf_score_attr, tmp_dir):
    # output files
    compare_file = os.path.join(output_dir, "compare_transcripts.txt")
    global_stats_file = os.path.join(output_dir, "global_stats.txt")
    tmp_gtf_file = os.path.join(output_dir, "tmp.gtf")
    tmp_sorted_gtf_file = os.path.splitext(tmp_gtf_file)[0] + ".srt.gtf"
    # merge and sort ref/test gtf files
    logging.info("Merging reference and test GTF files")
    # make temporary file to store merged ref/test gtf files
    outfh = open(tmp_gtf_file, "w")
    logging.info("Adding reference GTF file")
    add_gtf_file(ref_gtf_file, outfh, is_ref=True, sample_id=None)
    logging.info("Adding test GTF file")
    add_gtf_file(test_gtf_file, outfh, is_ref=False, sample_id='assembly')
    outfh.close()
    logging.info("Sorting merged GTF file")
    sort_gtf(tmp_gtf_file, tmp_sorted_gtf_file, tmp_dir=tmp_dir)
    os.remove(tmp_gtf_file)
    # compare assemblies
    logging.info("Comparing assemblies")
    cmp_fh = open(compare_file, "w")
    print >>cmp_fh, '\t'.join(map(str, MatchStats.header_fields()))
    stats_obj = GlobalStats()
    for locus_transcripts in parse_gtf(open(tmp_sorted_gtf_file)):
        locus_chrom = locus_transcripts[0].chrom
        locus_start = locus_transcripts[0].start
        locus_end = max(t.end for t in locus_transcripts)
        logging.debug("[LOCUS] %s:%d-%d %d transcripts" % 
                      (locus_chrom, locus_start, locus_end, 
                       len(locus_transcripts)))    
        # score transcripts
        for t in locus_transcripts:
            if gtf_score_attr is None:
                t.score = 0.0
            else:
                t.score = float(t.attrs.get(gtf_score_attr, 0.0))
        # run comparison
        for mobj in compare_locus(locus_transcripts):
            print >>cmp_fh, str(mobj)
        # measure global stats
        locus_stats_obj = gather_global_stats(locus_transcripts)
        stats_obj = stats_obj + locus_stats_obj        
    # cleanup
    cmp_fh.close()
    logging.info("Printing report")    
    f = open(global_stats_file, "w")
    print >>f, stats_obj.report()
    os.remove(tmp_sorted_gtf_file)
    logging.info("Done")    


def main():
    # parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("-v", "--verbose", action="store_true", 
                        dest="verbose", default=False)
    parser.add_argument("-o", "--output-dir", dest="output_dir", 
                        default="compare")
    parser.add_argument("--gtf-score-attr", dest="gtf_score_attr", 
                        default=None)
    parser.add_argument('--tmp-dir', dest="tmp_dir", default=None,
                        help="directory for sort to store temp files")
    parser.add_argument("ref_gtf_file")
    parser.add_argument("test_gtf_file")
    args = parser.parse_args()
    # check command line parameters
    if not os.path.exists(args.ref_gtf_file):
        parser.error("Reference GTF file %s not found" % (args.ref_gtf_file))
    if not os.path.exists(args.test_gtf_file):
        parser.error("Test GTF file %s not found" % (args.test_gtf_file))
    if not os.path.exists(args.output_dir):
        logging.debug("Creating output directory '%s'" % (args.output_dir))
        os.makedirs(args.output_dir)
    # set logging level
    if args.verbose:
        level = logging.DEBUG
    else:
        level = logging.INFO
    logging.basicConfig(level=level,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    logging.info("AssemblyLine %s" % (assemblyline.__version__))
    logging.info("----------------------------------")   
    # show parameters
    logging.info("Parameters:")
    logging.info("verbose logging:       %s" % (args.verbose))
    logging.info("reference gtf file:    %s" % (args.ref_gtf_file))
    logging.info("test gtf file:         %s" % (args.test_gtf_file))
    logging.info("output directory:      %s" % (args.output_dir))
    logging.info("test gtf score attr:   %s" % (args.gtf_score_attr))
    compare_assembly(args.ref_gtf_file, args.test_gtf_file,
                     args.output_dir, args.gtf_score_attr,
                     args.tmp_dir)
    return 0

if __name__ == "__main__":
    sys.exit(main())