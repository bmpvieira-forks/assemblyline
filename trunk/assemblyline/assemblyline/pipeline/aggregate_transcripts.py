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
from assemblyline.lib.librarytable import LibraryInfo
from assemblyline.lib.gtf import GTFFeature, sort_gtf
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
        feature.attrs[GTFAttr.COHORT_ID] = s.cohort_id
        feature.attrs[GTFAttr.SAMPLE_ID] = s.sample_id
        feature.attrs[GTFAttr.LIBRARY_ID] = s.library_id
        feature.attrs[GTFAttr.REF] = '0'
        if rename_ids:
            t_id = feature.attrs[GTFAttr.TRANSCRIPT_ID]
            g_id = feature.attrs[GTFAttr.GENE_ID]
            new_t_id = "%s.%s" % (s.library_id, t_id)
            new_g_id = "%s.%s" % (s.library_id, g_id)
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
    parser.add_argument('ref_gtf_file')
    parser.add_argument('sample_table_file')
    parser.add_argument('output_file')
    args = parser.parse_args()
    # check command line parameters
    if not os.path.exists(args.ref_gtf_file):
        parser.error("Reference GTF file %s not found" % (args.ref_gtf_file))
    if not os.path.exists(args.sample_table_file):
        parser.error("Sample table file %s not found" % (args.sample_table_file))
    # parse sample table
    logging.info("Parsing sample table")
    libinfos = []
    valid = True
    for lib in LibraryInfo.from_file(args.sample_table_file):
        # exclude samples
        if not lib.is_valid():
            logging.error("\tcohort=%s patient=%s sample=%s library=%s not valid" % 
                           (lib.cohort_id, lib.patient_id, lib.sample_id, lib.library_id))
            valid = False
        else:
            libinfos.append(lib)
    if not valid:
        logging.warning("Invalid samples in library table file")
    # show parameters
    logging.info("Parameters:")
    logging.info("Reference GTF file:   %s" % (args.ref_gtf_file))
    logging.info("Sample table file:    %s" % (args.sample_table_file))
    logging.info("Output file:          %s" % (args.output_file))
    # open tmp file
    tmp_file = os.path.splitext(args.output_file)[0] + ".unsorted.gtf"
    tmp_fh = open(tmp_file, "w")
    # read reference GTF file and bin by transcript id
    logging.info("Adding reference GTF file")
    add_reference_gtf_file(args.ref_gtf_file, tmp_fh)
    # parse sample table
    logging.info("Adding samples")
    for lib in libinfos:
        logging.info("\tadding cohort=%s patient=%s sample=%s library=%s" % 
                     (lib.cohort_id, lib.patient_id, lib.sample_id, lib.library_id))
        add_sample_gtf_file(lib, tmp_fh, args.rename_ids)
    tmp_fh.close()
    logging.info("Sorting GTF")
    sort_gtf(tmp_file, args.output_file)
    os.remove(tmp_file)
    return 0

if __name__ == '__main__':
    sys.exit(main())
