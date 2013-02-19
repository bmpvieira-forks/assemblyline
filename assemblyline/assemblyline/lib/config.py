'''
Created on Feb 11, 2013

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
import os

# file I/O parameters
MAX_OPEN_FILE_DESCRIPTORS = 1000 

# default parameter settings
MIN_TRANSCRIPT_LENGTH = 250
MIN_EXON_LENGTH = 15
CLIP_HIGH_PERCENTILE = 99.5

# prior observations of rna bases in genomic regions
PRIOR_MRNA_PROB = 0.811
PRIOR_INTRON_PROB = 0.1148
PRIOR_INTERGENIC_PROB = (PRIOR_MRNA_PROB - PRIOR_INTRON_PROB)
MODEL_ANN_COV_RATIO = 0.90
MODEL_ANN_INTRON_RATIO = 0.90

# percentiles for transcript score statistics
TRANSCRIPT_STATS_FIELDS = ['passed', 'failed', 'too_low_score', 
                           'too_short',  'too_short_exon', 
                           'clip_high_score', 'clipped'] 
TRANSCRIPT_SCORE_QUANTILES = [0, 1, 5, 10, 25, 50, 75, 90, 95, 99, 100]

# maximum bp for nearest gene distance calculation
MAX_GENE_DIST = 1000000
CATEGORY_PROTEIN = "protein"
CATEGORY_NCRNA = "ncrna"
CATEGORY_ANTISENSE = "antisense"
CATEGORY_INTRONIC = "intronic"
CATEGORY_INTERGENIC = "intergenic"

# modes for scoring input transcripts
SCORING_MODES = ("unweighted", "gtf_attr")

# constant file names
LIBRARY_ID_MAP = 'library_id.map'
SAMPLE_ID_MAP = 'sample_id.map'
TRANSCRIPTS_DROPPED_GTF_FILE = "transcripts.dropped.gtf"
TRANSCRIPTS_GTF_FILE = "transcripts.gtf"
TRANSCRIPT_STATS_FILE = "transcript_stats.txt"
ANNOTATED_TRANSCRIPTS_GTF_FILE = 'transcripts.annotated.gtf'
CATEGORY_COUNTS_FILE = "category_counts.txt"
CLASSIFY_DIR = 'classify'
REF_GTF_FILE = 'ref.gtf'
CLASSIFY_REPORT_FILE = "classification_report.txt"
BACKGROUND_GTF_FILE = 'background.gtf'
EXPRESSED_GTF_FILE = 'expressed.gtf'

class AssemblylineResults(object):
    def __init__(self, output_dir):
        self.run_dir = output_dir
        self.tmp_dir = os.path.join(output_dir, 'tmp')
        self.library_id_map = os.path.join(output_dir, LIBRARY_ID_MAP)
        self.sample_id_map = os.path.join(output_dir, SAMPLE_ID_MAP)
        self.transcripts_dropped_gtf_file = os.path.join(output_dir, TRANSCRIPTS_DROPPED_GTF_FILE)
        self.transcripts_gtf_file = os.path.join(output_dir, TRANSCRIPTS_GTF_FILE)
        self.transcript_stats_file = os.path.join(output_dir, TRANSCRIPT_STATS_FILE)
        self.annotated_transcripts_gtf_file = os.path.join(output_dir, ANNOTATED_TRANSCRIPTS_GTF_FILE)
        self.classify_dir = os.path.join(output_dir, CLASSIFY_DIR)
        self.category_counts_file = os.path.join(output_dir, CATEGORY_COUNTS_FILE)
        self.ref_gtf_file = os.path.join(output_dir, REF_GTF_FILE)
        self.background_gtf_file = os.path.join(output_dir, BACKGROUND_GTF_FILE)
        self.expressed_gtf_file = os.path.join(output_dir, EXPRESSED_GTF_FILE)
        self.classify_report_file = os.path.join(output_dir, CLASSIFY_REPORT_FILE)
        