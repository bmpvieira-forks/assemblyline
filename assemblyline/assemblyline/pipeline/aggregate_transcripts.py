'''
Created on Nov 28, 2011

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
import logging
import argparse
import sys
import os
import collections
import operator

import assemblyline
from assemblyline.lib.sampletable import SampleInfo
from assemblyline.lib.gtf import GTFFeature
from assemblyline.lib.base import GTFAttr

def add_reference_gtf_file(ref_gtf_file, outfh):
    transcript_dict = collections.defaultdict(lambda: [])
    for feature in GTFFeature.parse(open(ref_gtf_file)):
        if feature.feature_type != "exon":
            continue
        t_id = feature.attrs[GTFAttr.TRANSCRIPT_ID]
        transcript_dict[t_id].append(feature)
    # output reference transcripts
    for t_id, features in transcript_dict.iteritems():
        # sort features (exons) by start position
        features.sort(key=operator.attrgetter('start'))
        # annotate exons as reference features
        for f in features:
            f.attrs[GTFAttr.REF] = '1'
            print >>outfh, str(f)
        # transcript feature
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
        f.attrs[GTFAttr.REF] = '1'
        print >>outfh, str(f)

def add_sample_gtf_file(s, outfh, rename_ids):
    # parse and annotate GTF files
    for feature in GTFFeature.parse(open(s.gtf_file)):
        feature.attrs[GTFAttr.COHORT_ID] = s.cohort
        feature.attrs[GTFAttr.SAMPLE_ID] = s.sample
        feature.attrs[GTFAttr.LIBRARY_ID] = s.library
        feature.attrs[GTFAttr.REF] = '0'
        if rename_ids:
            t_id = feature.attrs[GTFAttr.TRANSCRIPT_ID]
            g_id = feature.attrs[GTFAttr.GENE_ID]
            new_t_id = "%s.%s" % (s.library, t_id)
            new_g_id = "%s.%s" % (s.library, g_id)
            feature.attrs[GTFAttr.TRANSCRIPT_ID] = new_t_id
            feature.attrs[GTFAttr.GENE_ID] = new_g_id
        print >>outfh, str(feature)

def main():
    # setup logging
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    logging.info("AssemblyLine %s" % (assemblyline.__version__))
    logging.info("----------------------------------")
    # parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument('--rename-ids', dest="rename_ids", 
                        action="store_true", default=False,
                        help="prepend library id to transcript and gene id "
                        "attributes to avoid mistaking transcripts from "
                        "different libraries as the same because of "
                        "redundant transcript ids")
    parser.add_argument('-o', dest="output_file", default=None)
    parser.add_argument('ref_gtf_file')
    parser.add_argument('sample_table_file')
    args = parser.parse_args()
    # check command line parameters
    if not os.path.exists(args.ref_gtf_file):
        parser.error("Reference GTF file %s not found" % (args.ref_gtf_file))
    if not os.path.exists(args.sample_table_file):
        parser.error("Sample table file %s not found" % (args.sample_table_file))
    # parse sample table
    logging.info("Parsing sample table")
    sampleinfos = []
    valid = True
    for s in SampleInfo.from_file(args.sample_table_file):
        # exclude samples
        if not s.is_valid():
            logging.error("\tcohort=%s patient=%s sample=%s library=%s not valid" % 
                          (s.cohort, s.patient, s.sample, s.library))
            valid = False
        else:
            sampleinfos.append(s)
    if not valid:
        parser.error("Invalid samples in sample table file")
    # show parameters
    logging.info("Reference GTF file:   %s" % (args.ref_gtf_file))
    logging.info("Sample table file:    %s" % (args.sample_table_file))
    if args.output_file is None:
        logging.info("Output file:          <stdout>")
    else:
        logging.info("Output file:          %s" % (args.output_file))
    # open output file
    if args.output_file is not None:
        outfh = open(args.output_file, "w")
    else:
        outfh = sys.stdout
    # read reference GTF file and bin by transcript id
    logging.info("Adding reference GTF file")
    add_reference_gtf_file(args.ref_gtf_file, outfh)
    # parse sample table
    logging.info("Adding samples")
    for s in sampleinfos:
        logging.info("\tadding cohort=%s patient=%s sample=%s library=%s" % 
                     (s.cohort, s.patient, s.sample, s.library))
        add_sample_gtf_file(s, outfh, args.rename_ids)
    # cleanup
    if args.output_file is not None:
        outfh.close()
    return 0

if __name__ == '__main__':
    sys.exit(main())
