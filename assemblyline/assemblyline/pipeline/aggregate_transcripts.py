'''
Created on Nov 28, 2011

@author: mkiyer

AssemblyLine: transcriptome meta-assembly from RNA-Seq

Copyright (C) 2012-2013 Matthew Iyer

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
import sys
import os
import collections
import operator
import random

# project imports
import assemblyline
import assemblyline.lib.config as config
from assemblyline.lib.base import Library, GTFAttr
from assemblyline.lib.gtf import GTFFeature, sort_gtf
from assemblyline.lib.stats import ECDF, scoreatpercentile

def make_transcript_feature(exon_features):
    f = GTFFeature()
    f.seqid = exon_features[0].seqid
    f.source = exon_features[0].source
    f.feature_type = 'transcript'
    f.start = exon_features[0].start
    f.end = exon_features[-1].end
    f.score = exon_features[0].score
    f.strand = exon_features[0].strand
    f.phase = '.'
    f.attrs = exon_features[0].attrs.copy()
    if "exon_number" in f.attrs:
        del f.attrs["exon_number"]
    return f

def add_reference_gtf_file(ref_gtf_file, test_gene_ids, 
                           random_test_frac, outfh):
    gene_dict = collections.defaultdict(lambda: [])
    user_defined_tests = len(test_gene_ids) > 0
    # group by gene id
    for feature in GTFFeature.parse(open(ref_gtf_file)):
        if feature.feature_type != "exon":
            continue
        # group by gene id        
        g_id = feature.attrs[GTFAttr.GENE_ID]
        gene_dict[g_id].append(feature)
    # output reference transcripts
    for g_id, g_features in gene_dict.iteritems():
        # label test transcripts
        if user_defined_tests:
            is_test = (g_id in test_gene_ids)
        else:
            is_test = (random.random() < random_test_frac)
        # group by transcript id
        transcript_dict = collections.defaultdict(lambda: [])
        for feature in g_features:
            t_id = feature.attrs[GTFAttr.TRANSCRIPT_ID]
            transcript_dict[t_id].append(feature)
        for t_id, t_features in transcript_dict.iteritems():
            # sort features (exons) by start position
            t_features.sort(key=operator.attrgetter('start'))
            # annotate exons as reference features
            for f in t_features:
                f.attrs[GTFAttr.REF] = '1'
                f.attrs[GTFAttr.TEST] = '1' if is_test else '0'
                print >>outfh, str(f)
            f = make_transcript_feature(t_features)
            f.attrs[GTFAttr.REF] = '1'
            f.attrs[GTFAttr.TEST] = '1' if is_test else '0'
            print >>outfh, str(f)
        del transcript_dict
    del gene_dict

def read_gtf_file(library, gtf_score_attr):
    # read all transcripts
    t_dict = collections.OrderedDict()
    cur_t_id = 1
    cur_g_id = 1
    t_id_map = {}
    g_id_map = {}
    for feature in GTFFeature.parse(open(library.gtf_file)):
        if feature.feature_type == "exon":
            t_id = feature.attrs[GTFAttr.TRANSCRIPT_ID]
            # rename transcript id
            if t_id not in t_id_map:
                new_t_id = "%s.T%d" % (library.library_id, cur_t_id)
                t_id_map[t_id] = new_t_id
                cur_t_id += 1
            else:
                new_t_id = t_id_map[t_id]
            # rename gene id
            g_id = feature.attrs[GTFAttr.GENE_ID]
            if g_id not in g_id_map:
                new_g_id = "%s.G%d" % (library.library_id, cur_g_id)
                g_id_map[g_id] = new_g_id
                cur_g_id += 1
            else:
                new_g_id = g_id_map[g_id]
            # update transcript attributes
            newattrs = {GTFAttr.TRANSCRIPT_ID: new_t_id,
                        GTFAttr.GENE_ID: new_g_id,
                        GTFAttr.SAMPLE_ID: library.sample_id,
                        GTFAttr.LIBRARY_ID: library.library_id,
                        GTFAttr.REF: '0',
                        GTFAttr.SCORE: feature.attrs.get(gtf_score_attr, '0.0')}
            feature.attrs = newattrs
            # store feature
            if new_t_id not in t_dict:
                t_dict[new_t_id] = []
            t_dict[new_t_id].append(feature)
    return t_dict

def filter_transcripts(library_id, t_dict, outfileh, dropfileh, statsfileh, 
                       min_transcript_length):
    # filter transcripts
    passed = 0
    failed = 0
    too_short = 0
    too_short_exon = 0
    passed_scores = []
    failed_scores = []
    filtered_t_dict = collections.OrderedDict()
    for t_id, features in t_dict.iteritems():
        strand = features[0].strand
        reverse = True if strand == "-" else False
        features.sort(key=operator.attrgetter('start'), reverse=reverse)
        # check first/last exon lengths and clip very short exons
        num_exons = len(features)
        if num_exons > 1:
            f = features[0]
            length = f.end - f.start
            if length < config.MIN_EXON_LENGTH:
                too_short_exon += 1
                features = features[1:]
            f = features[-1]
            length = f.end - f.start
            if length < config.MIN_EXON_LENGTH:
                too_short_exon += 1
                features = features[:-1]
        # filter transcripts by various criteria
        score = float(features[0].attrs[GTFAttr.SCORE])
        transcript_length = sum((f.end - f.start) for f in features)
        keep = True
        if (transcript_length <= min_transcript_length):
            too_short += 1
            keep = False
        if not keep:
            failed += 1
            failed_scores.append(score)
            feature = make_transcript_feature(features)
            print >>dropfileh, str(feature)
            for feature in features:
                print >>dropfileh, str(feature)
        else:
            passed += 1
            passed_scores.append(score)
            filtered_t_dict[t_id] = (score, features) 
    # make empirical cdf for scores
    ecdf = ECDF(passed_scores, side="left")
    for t_id, vtuple in filtered_t_dict.iteritems():
        score, features = vtuple
        pctrank = 100.0 * ecdf(score)
        # reverse features if this is negative strand
        if features[0].strand == "-":
            features.reverse()
        # write transcript
        feature = make_transcript_feature(features)
        feature.attrs[GTFAttr.SCORE] = score
        feature.attrs[GTFAttr.PCTRANK] = pctrank
        print >>outfileh, str(feature)
        # write exons
        for i,feature in enumerate(features):
            feature.attrs[GTFAttr.SCORE] = score
            feature.attrs[GTFAttr.PCTRANK] = pctrank
            feature.attrs['exon_number'] = "%d" % (i + 1)
            print >>outfileh, str(feature)

    # compute and write stats
    failed_quantiles = [scoreatpercentile(failed_scores, q) 
                        for q in config.TRANSCRIPT_SCORE_QUANTILES]
    passed_quantiles = [scoreatpercentile(passed_scores, q) 
                        for q in config.TRANSCRIPT_SCORE_QUANTILES]
    fields = [library_id, passed, failed, too_short, too_short_exon]
    fields.extend(failed_quantiles)
    fields.extend(passed_quantiles)
    print >>statsfileh, '\t'.join(map(str, fields))

def main():
    # setup logging
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    logging.info("AssemblyLine %s" % (assemblyline.__version__))
    logging.info("----------------------------------")
    # parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument('--min-transcript-length', type=int, 
                        dest="min_transcript_length",
                        metavar="N",
                        default=config.MIN_TRANSCRIPT_LENGTH,
                        help="Skip ab initio transcripts equal to or below "
                        "this length [default=%(default)s]")
    parser.add_argument("--gtf-score-attr", dest="gtf_score_attr", 
                        default="FPKM", metavar="ATTR",
                        help="GTF attribute field containing transcript "
                        "expression [default='%(default)s']")
    parser.add_argument('-o', '--output-dir', dest="output_dir", 
                        default="assemblyline_out",
                        help="directory to store assemblyline results and "
                        "intermediate files [default=%(default)s]")
    parser.add_argument("--random-test-frac", dest="random_test_frac", 
                        default=0.1, metavar="FRAC", type=float,
                        help="if no user-defined tests are specified "
                        "using '--tests' randomly designate a fraction "
                        "of reference transcripts as test data for use "
                        "in classification [default=%(default)s]")
    parser.add_argument("--tests", dest="test_file", default=None,
                        help="(optional) text file containing "
                        "reference 'gene_id' attributes "
                        "(one per line) that define test cases "
                        "to use for validation purposes")
    parser.add_argument('ref_gtf_file')
    parser.add_argument('library_table_file')
    args = parser.parse_args()
    # check command line parameters
    if not os.path.exists(args.library_table_file):
        parser.error("library table file %s not found" % (args.library_table_file))
    if args.min_transcript_length < 0:
        parser.error("min_transcript_length < 0")
    if not os.path.exists(args.ref_gtf_file):
        parser.error("reference GTF file %s not found" % (args.ref_gtf_file))
    if (args.test_file is not None) and (not os.path.exists(args.test_file)):
        parser.error("test file %s not found" % (args.test_file))
    if (args.random_test_frac < 0):
        parser.error("cannot set --random-test-frac < 0")
    # show parameters
    logging.info("Parameters:")
    logging.info("min transcript length: %d" % (args.min_transcript_length))
    logging.info("gtf score attr:        %s" % (args.gtf_score_attr))
    logging.info("output directory:      %s" % (args.output_dir))
    logging.info("reference GTF file:    %s" % (args.ref_gtf_file))
    logging.info("test file:             %s" % (args.test_file))
    logging.info("library table file:    %s" % (args.library_table_file))
    logging.info("----------------------------------")
    # setup results
    results = config.AssemblylineResults(args.output_dir)
    # create output directory
    if not os.path.exists(results.run_dir):
        logging.debug("Creating output directory '%s'" % (results.run_dir))
        os.makedirs(results.run_dir)
    if not os.path.exists(results.tmp_dir):
        logging.info("Creating tmp directory '%s'" % (results.tmp_dir))
        os.makedirs(results.tmp_dir)
    # parse sample table
    logging.info("Parsing library table")
    libraries = []
    valid = True
    library_num = 1
    sample_num = 1
    sample_id_map = {}
    library_map_fileh = open(results.library_id_map, 'w')
    sample_map_fileh = open(results.sample_id_map, 'w')
    for library in Library.from_file(args.library_table_file):
        # exclude samples
        if not os.path.exists(library.gtf_file):
            logging.warning("Library '%s' GTF file not found" % (library.library_id)) 
            continue
        # rename library id
        new_library_id = "L%d" % (library_num)
        print >>library_map_fileh, '\t'.join([new_library_id, library.library_id]) 
        library_num += 1
        library.library_id = new_library_id
        # rename sample id
        if library.sample_id not in sample_id_map:
            new_sample_id = "S%d" % (sample_num)
            print >>sample_map_fileh, '\t'.join([new_sample_id, library.sample_id]) 
            sample_id_map[library.sample_id] = new_sample_id
            sample_num += 1
        else:
            new_sample_id = sample_id_map[library.sample_id]
        library.sample_id = new_sample_id
        libraries.append(library)
    if not valid:
        logging.warning("Invalid libraries in table file")
    library_map_fileh.close()
    sample_map_fileh.close()
    # setup output files
    tmp_file = os.path.join(results.tmp_dir, "transcripts.unsorted.gtf")
    tmpfileh = open(tmp_file, "w")
    dropfileh = open(results.transcripts_dropped_gtf_file, "w")
    statsfileh = open(results.transcript_stats_file, 'w')
    header_fields = ['#library_id']
    header_fields.extend(config.TRANSCRIPT_STATS_FIELDS)
    header_fields.extend([("failed_q%d" % x) for x in config.TRANSCRIPT_SCORE_QUANTILES])
    header_fields.extend([("passed_q%d" % x) for x in config.TRANSCRIPT_SCORE_QUANTILES])
    print >>statsfileh, '\t'.join(header_fields)
    # read test transcripts
    test_gene_ids = set()
    if args.test_file is not None:
        fileh = open(args.test_file)
        test_gene_ids.update(line.strip() for line in fileh)
        fileh.close()
        logging.info("Read %d test genes" % len(test_gene_ids))
    # read reference GTF file and aggregate
    logging.info("Adding reference GTF file")
    add_reference_gtf_file(args.ref_gtf_file, test_gene_ids, 
                           args.random_test_frac, tmpfileh)
    # parse sample table
    logging.info("Adding libraries")
    for library in libraries:
        t_dict = read_gtf_file(library, args.gtf_score_attr)
        logging.debug("Read %s transcripts from file %s" % (len(t_dict), 
                                                            library.gtf_file))
        if len(t_dict) == 0:
            logging.warning("Library %s has no transcripts" % 
                            (library.library_id))
        else:
            filter_transcripts(library.library_id, t_dict, 
                               tmpfileh, dropfileh, statsfileh, 
                               args.min_transcript_length)
    statsfileh.close()
    tmpfileh.close()
    logging.info("Sorting GTF")
    retcode = sort_gtf(tmp_file, results.transcripts_gtf_file, tmp_dir=results.tmp_dir)
    if retcode != 0:
        logging.error("sort GTF failed")
        if os.path.exists(results.transcripts_gtf_file):
            os.remove(results.transcripts_gtf_file)
    os.remove(tmp_file)
    logging.info("Done")
    return retcode

if __name__ == '__main__':
    sys.exit(main())
