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
import tempfile
import operator
import collections

import assemblyline
from assemblyline.lib.gtf import GTFFeature, sort_gtf
from assemblyline.lib.base import GTFAttr
from assemblyline.lib.transcript_parser import parse_gtf
from assemblyline.lib.transcript import NO_STRAND, POS_STRAND, NEG_STRAND
from assemblyline.lib.assemble.transcript_graph import create_undirected_transcript_graph, get_transcript_node_map
from assemblyline.lib.assemble.base import STRAND_SCORE, IS_REF
from assemblyline.pipeline.aggregate_transcripts import add_reference_gtf_file

class CompareData(object):
    __slots__ = ('is_ref', 'is_assembly')
    def __init__(self):
        self.is_ref = False
        self.is_assembly = False

class CompareStats(object):
    fields = ('introns_both', 'introns_ref_only', 'introns_assembly_only',
              'patterns_both', 'patterns_ref_only', 'patterns_assembly_only', 
              'cov_both', 'cov_ref_only', 'cov_assembly_only')
    def __init__(self):
        self.introns_both = 0
        self.introns_ref_only = 0
        self.introns_assembly_only = 0
        self.patterns_both = 0
        self.patterns_ref_only = 0
        self.patterns_assembly_only = 0
        self.cov_both = 0
        self.cov_ref_only = 0
        self.cov_assembly_only = 0
    def __add__(self, x):
        z = CompareStats()
        for f in CompareStats.fields:
            val = getattr(self, f) + getattr(x, f)
            setattr(z, f, val)
        return z

class TranscriptStats(object):
    def __init__(self):
        self.transcript_id = None
        self.gene_id = None
        self.score = None
        self.length = None
        self.found_length = None
        self.num_introns = None
        self.found_introns = None
        self.found_pattern = None
    @staticmethod
    def header_fields():
        return ["gene_id", "transcript_id", "score", "length",
                "assembly_length", "num_introns", "assembly_introns",
                "assembly_pattern"]        
    def __str__(self):
        fields = [self.gene_id, self.transcript_id, self.score, self.length,
                  self.found_length, self.num_introns, self.found_introns,
                  int(self.found_pattern)]
        return '\t'.join(map(str,fields))
    
def measure_transcript_sensitivity(ref_transcripts, test_transcripts):
    strand_introns = {POS_STRAND: set(),
                      NEG_STRAND: set()}
    strand_patterns = {POS_STRAND: set(),
                       NEG_STRAND: set()}
    for t in test_transcripts:
        if t.strand == NO_STRAND:
            continue
        pattern = []
        for start,end in t.iterintrons():
            strand_introns[t.strand].add((start,end))
            pattern.append((start,end))
        strand_patterns[t.strand].add(tuple(pattern))
    # add all transcripts to a graph
    all_transcripts = []
    all_transcripts.extend(ref_transcripts)
    all_transcripts.extend(test_transcripts)
    G = create_undirected_transcript_graph(all_transcripts)
    # build a mapping from transcripts to graph nodes using the 
    # transcript id attributes of the nodes
    transcript_node_map = get_transcript_node_map(G)
    # compile statistics
    for t in ref_transcripts:
        found_introns = 0
        total_introns = 0
        pattern = []
        for start,end in t.iterintrons():
            if (start,end) in strand_introns[t.strand]:
                found_introns += 1
            total_introns += 1
            pattern.append((start,end))
        found_pattern = tuple(pattern) in strand_patterns[t.strand]
        total_score = 0.0
        total_length = 0
        found_length = 0
        nodes = transcript_node_map[t.attrs[GTFAttr.TRANSCRIPT_ID]]
        for n in nodes:
            length = (n.end - n.start)
            score = G.node[n][STRAND_SCORE][t.strand]
            if score > 0.0:
                found_length += length
            total_score += score
            total_length += length
        ts = TranscriptStats()
        ts.transcript_id = t.attrs[GTFAttr.TRANSCRIPT_ID]
        ts.gene_id = t.attrs[GTFAttr.GENE_ID]
        ts.score = total_score / float(total_length)
        ts.length = total_length
        ts.found_length = found_length
        ts.num_introns = total_introns
        ts.found_introns = found_introns
        ts.found_pattern = found_pattern
        yield ts

def compare_locus(transcripts, gtf_score_attr):
    # handle splicing related statistics
    strand_intron_dict = {POS_STRAND: collections.defaultdict(lambda: CompareData()),
                          NEG_STRAND: collections.defaultdict(lambda: CompareData())}
    strand_splicing_pattern_dict = {POS_STRAND: collections.defaultdict(lambda: CompareData()),
                                    NEG_STRAND: collections.defaultdict(lambda: CompareData())}
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
        d = pattern_dict[tuple(pattern)]
        if is_ref:
            d.is_ref = True
        else:
            d.is_assembly = True
    # compile statistics
    stats_obj = CompareStats()
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
    G = create_undirected_transcript_graph(transcripts)
    for n,nd in G.nodes_iter(data=True):
        length = (n.end - n.start)
        strand_is_ref = nd[IS_REF]
        strand_score = nd[STRAND_SCORE]
        for strand in xrange(0,3):
            is_ref = strand_is_ref[strand]
            score = strand_score[strand]
            is_assembly = (score > 0.0)
            if is_ref and is_assembly:
                stats_obj.cov_both += length
            elif is_ref:
                stats_obj.cov_ref_only += length
            elif is_assembly:
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
                     gtf_score_attr):
    # output files
    ref_output_file = os.path.join(output_dir, "ref_transcripts.txt")
    stats_file = os.path.join(output_dir, "stats.txt")
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
    sort_gtf(tmp_gtf_file, tmp_sorted_gtf_file)
    os.remove(tmp_gtf_file)
    # start comparing assemblies
    logging.info("Comparing assemblies")
    ref_fh = open(ref_output_file, "w")
    print >>ref_fh, '\t'.join(map(str, TranscriptStats.header_fields()))
    stats_obj = CompareStats()
    for locus_transcripts in parse_gtf(open(tmp_sorted_gtf_file)):
        locus_chrom = locus_transcripts[0].chrom
        locus_start = locus_transcripts[0].start
        locus_end = max(t.end for t in locus_transcripts)
        logging.debug("[LOCUS] %s:%d-%d %d transcripts" % 
                      (locus_chrom, locus_start, locus_end, 
                       len(locus_transcripts)))
        # score and divide into ref/test transcript lists
        ref_transcripts = []
        test_transcripts = []
        for t in locus_transcripts:
            is_ref = bool(int(t.attrs[GTFAttr.REF]))
            if is_ref:
                t.score = 0.0
                ref_transcripts.append(t)
            else:
                t.score = float(t.attrs[gtf_score_attr])
                test_transcripts.append(t)
        # measure sensitivity for detecting reference transcripts
        for ts in measure_transcript_sensitivity(ref_transcripts, 
                                                 test_transcripts):
            print >>ref_fh, str(ts)
        # measure overall stats
        locus_stats_obj = compare_locus(locus_transcripts, gtf_score_attr)
        stats_obj = stats_obj + locus_stats_obj
    ref_fh.close()
    # print stats report
    introns_total = stats_obj.introns_both + stats_obj.introns_ref_only + stats_obj.introns_assembly_only
    patterns_total = stats_obj.patterns_both + stats_obj.patterns_ref_only + stats_obj.patterns_assembly_only
    cov_total = stats_obj.cov_both + stats_obj.cov_ref_only + stats_obj.cov_assembly_only
    f = open(stats_file, "w")
    print >>f, "introns_total=%d" % (introns_total)
    print >>f, "introns_both=%d" % (stats_obj.introns_both)
    print >>f, "introns_ref_only=%d" % (stats_obj.introns_ref_only)
    print >>f, "introns_assembly_only=%d" % (stats_obj.introns_assembly_only)
    print >>f, "introns_precision=%f" % (stats_obj.introns_both / float(stats_obj.introns_both + stats_obj.introns_assembly_only))
    print >>f, "introns_recall=%f" % (stats_obj.introns_both / float(stats_obj.introns_both + stats_obj.introns_ref_only))
    print >>f, "patterns_total=%d" % (patterns_total)
    print >>f, "patterns_both=%d" % (stats_obj.patterns_both)
    print >>f, "patterns_ref_only=%d" % (stats_obj.patterns_ref_only)
    print >>f, "patterns_assembly_only=%d" % (stats_obj.patterns_assembly_only)
    print >>f, "patterns_precision=%f" % (stats_obj.patterns_both / float(stats_obj.patterns_both + stats_obj.patterns_assembly_only))
    print >>f, "patterns_recall=%f" % (stats_obj.patterns_both / float(stats_obj.patterns_both + stats_obj.patterns_ref_only))
    print >>f, "cov_total=%d" % (cov_total)
    print >>f, "cov_both=%d" % (stats_obj.cov_both)
    print >>f, "cov_ref_only=%d" % (stats_obj.cov_ref_only)
    print >>f, "cov_assembly_only=%d" % (stats_obj.cov_assembly_only)
    print >>f, "cov_precision=%f" % (stats_obj.cov_both / float(stats_obj.cov_both + stats_obj.cov_assembly_only))
    print >>f, "cov_recall=%f" % (stats_obj.cov_both / float(stats_obj.cov_both + stats_obj.cov_ref_only))
    f.close()
    # cleanup
    os.remove(tmp_sorted_gtf_file)

def main():
    # parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("-v", "--verbose", action="store_true", 
                        dest="verbose", default=False)
    parser.add_argument("-o", "--output-dir", dest="output_dir", 
                        default="compare")
    parser.add_argument("--gtf-score-attr", dest="gtf_score_attr", 
                        default="score")
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
                     args.output_dir, args.gtf_score_attr)
    return 0

if __name__ == "__main__":
    sys.exit(main())