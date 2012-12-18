'''
Created on Aug 3, 2011

@author: mkiyer
'''
import re
import os
import logging
import collections
import xml.etree.cElementTree as etree

from assemblyline.rnaseq.lib.base import check_executable, check_sam_file, indent_xml, file_exists_and_nz_size, parse_bool
from assemblyline.rnaseq.lib.libtable import FRAGMENT_LAYOUT_PAIRED, FR_UNSTRANDED
from assemblyline.rnaseq.lib.inspect import RnaseqLibraryMetrics

# default parameter values
MIN_INSPECT_SAMPLES = 100
MAX_INSPECT_SAMPLES = 1000000
STRAND_SPECIFIC_CUTOFF_FRAC = 0.90
# job return codes
JOB_SUCCESS = 0
JOB_ERROR = 1
# gzip/bz2 compress level
DEFAULT_COMPRESS_LEVEL = 5
# XML files
LIBRARY_XML_FILE = "library.xml"
CONFIG_XML_FILE = "config.xml"
# fastq file names
FASTQ_PREFIXES = ("read1", "read2")
FASTQ_FILES = [("%s.fq" % x) for x in FASTQ_PREFIXES] 
# fastqc
FASTQC_DIR_EXTENSION = "_fastqc"
FASTQC_DATA_FILE = "fastqc_data.txt"
FASTQC_REPORT_FILE = "fastqc_report.html"
# abundant sequence mapping
ABUNDANT_SAM_FILES = ('abundant_hits_read1.sam', 'abundant_hits_read2.sam')
ABUNDANT_BAM_FILE = 'abundant_hits.bam'
SORTED_ABUNDANT_BAM_FILE = 'abundant_hits.srt.bam'
ABUNDANT_COUNTS_FILE = 'abundant_counts.txt'
# filtered fastq files
FILTERED_FASTQ_PREFIX = 'filtered_read'
FILTERED_FASTQ_FILES = tuple(("%s%d.fq" % (FILTERED_FASTQ_PREFIX,x)) for x in (1,2)) 
# rnaseq library inspection
LIBRARY_METRICS_FILE = "library_metrics.txt"
FRAG_SIZE_DIST_PLOT_FILE = "frag_size_dist_plot.pdf"
# repeat elements
REPEAT_ELEMENTS_FILE = "repeat_element_counts.txt"
# tophat alignment results
TOPHAT_DIR = 'tophat'
TOPHAT_BAM_FILE = os.path.join(TOPHAT_DIR, "accepted_hits.bam")
TOPHAT_BAM_INDEX_FILE = TOPHAT_BAM_FILE + ".bai"
TOPHAT_JUNCTIONS_FILE = os.path.join(TOPHAT_DIR, "junctions.bed")
# tophat fusion results
TOPHAT_FUSION_DIR = 'tophatfusion'
TOPHAT_FUSION_BAM_FILE = os.path.join(TOPHAT_FUSION_DIR, "accepted_hits.bam")
TOPHAT_FUSION_BAM_INDEX_FILE = TOPHAT_FUSION_BAM_FILE + ".bai"
TOPHAT_FUSION_FILE = os.path.join(TOPHAT_FUSION_DIR, "fusions.out")
TOPHAT_FUSION_POST_RESULT_FILE = os.path.join(TOPHAT_FUSION_DIR, 'result.txt')
TOPHAT_FUSION_POST_HTML_FILE = os.path.join(TOPHAT_FUSION_DIR, 'result.html')
TOPHAT_FUSION_TMP_FILES = (os.path.join(TOPHAT_FUSION_DIR, f) 
                           for f in ("check", "fusion_seq.bwtout", 
                                     "fusion_seq.map", 
                                     "fusion_seq.fa", 
                                     "potential_fusion.txt"))
# picard metrics files
PICARD_ALIGNMENT_SUMMARY_METRICS = "picard.alignment_summary_metrics"
PICARD_INSERT_SIZE_HISTOGRAM_PDF = "picard.insert_size_histogram.pdf"
PICARD_INSERT_SIZE_METRICS = "picard.insert_size_metrics"
PICARD_QUALITY_BY_CYCLE_METRICS = "picard.quality_by_cycle_metrics"
PICARD_QUALITY_BY_CYCLE_PDF = "picard.quality_by_cycle.pdf"
PICARD_QUALITY_DISTRIBUTION_METRICS = "picard.quality_distribution_metrics"
PICARD_QUALITY_DISTRIBUTION_PDF = "picard.quality_distribution.pdf"
PICARD_RNASEQ_METRICS = "picard.rnaseq_metrics"
PICARD_RNASEQ_METRICS_PLOT_PDF = "picard.rnaseq_metrics_plot.pdf"
# coverage bedgraph file
COVERAGE_BIGWIG_PREFIX = "coverage"
COVERAGE_UCSC_TRACK_FILE = "coverage.ucsc.tracks"
# cufflinks output
CUFFLINKS_AB_INITIO_DIR = "cufflinks_ab_initio"
CUFFLINKS_KNOWN_DIR = "cufflinks_known"
CUFFLINKS_TRANSCRIPTS_GTF_FILE = "transcripts.gtf"
CUFFLINKS_GENES_FILE = "genes.fpkm_tracking"
CUFFLINKS_ISOFORMS_FILE = "isoforms.fpkm_tracking"
# htseq count output
HTSEQ_COUNT_KNOWN_OUTPUT_FILE = "htseq_count_known_genes.txt"
# picard mark duplicates output
TOPHAT_RMDUP_BAM_FILE = "accepted_hits.rmdup.bam"
PICARD_DUPLICATE_METRICS = "picard.duplicate_metrics"
# varscan output
VARSCAN_SNV_FILE = "varscan_snvs.vcf"
# job complete
JOB_DONE_FILE = "job.done"
# job memory and runtime
PBS_JOB_MEM = 24000
PBS_JOB_WALLTIME = "160:00:00"
# quality score formats
SANGER_FORMAT = "sanger"
SOLEXA_FORMAT = "solexa"
ILLUMINA_FORMAT = "illumina"
FASTQ_QUAL_FORMATS = (SANGER_FORMAT, SOLEXA_FORMAT, ILLUMINA_FORMAT)

class PipelineConfigError(Exception):
    pass

def check_tophat_juncs_file(filename):
    if not file_exists_and_nz_size(filename):
        return False
    lines = 0
    for line in open(filename):
        if not line or line.startswith("track"):
            continue
        lines += 1
        if lines > 100:
            return True
    return lines > 0

def check_library_metrics(filename):
    is_valid=True
    if not file_exists_and_nz_size(filename):
        is_valid = False
    else:
        try:
            obj = RnaseqLibraryMetrics.from_file(open(filename))
        except:
            is_valid = False    
    return is_valid

class RnaseqResults(object):
    def __init__(self, library, output_dir):
        self.library_id = library.library_id
        self.output_dir = output_dir
        self.tmp_dir = os.path.join(self.output_dir, "tmp")
        self.log_dir = os.path.join(self.output_dir, "log")
        # pipeline config used to run
        self.library_xml_file = os.path.join(self.output_dir, LIBRARY_XML_FILE)
        self.config_xml_file = os.path.join(self.output_dir, CONFIG_XML_FILE)
        # BAM FASTQ files
        self.bam_fastq_prefixes = []
        self.bam_read1_files = []
        self.bam_read2_files = []
        for i in xrange(1,len(library.bam_files)+1):
            prefix = os.path.join(self.tmp_dir, "bam_file%03d" % (i))
            self.bam_fastq_prefixes.append(prefix)
            self.bam_read1_files.append(prefix + "_1.fq")
            if library.fragment_layout == FRAGMENT_LAYOUT_PAIRED:
                self.bam_read2_files.append(prefix + "_2.fq")
        # FASTQ files
        self.copied_fastq_files = []
        if library.fragment_layout == FRAGMENT_LAYOUT_PAIRED:
            self.copied_fastq_files = [os.path.join(self.tmp_dir, x) for x in FASTQ_FILES]
        else:
            self.copied_fastq_files = [os.path.join(self.tmp_dir, FASTQ_FILES[0])]
        # FASTQC results
        self.fastqc_data_files = []
        self.fastqc_report_files = []
        for readnum in xrange(len(self.copied_fastq_files)):
            fastqc_dir = os.path.join(self.output_dir, "%s%s" % (FASTQ_FILES[readnum], FASTQC_DIR_EXTENSION))
            self.fastqc_data_files.append(os.path.join(fastqc_dir, FASTQC_DATA_FILE))
            self.fastqc_report_files.append(os.path.join(fastqc_dir, FASTQC_REPORT_FILE))
        # Filtered abundant BAM and FASTQ
        self.filtered_fastq_files = []
        for readnum in xrange(len(self.copied_fastq_files)):
            self.filtered_fastq_files.append(os.path.join(self.tmp_dir, FILTERED_FASTQ_FILES[readnum]))
        # sorted abundant reads bam file
        self.sorted_abundant_bam_file = os.path.join(self.output_dir, SORTED_ABUNDANT_BAM_FILE)
        self.abundant_counts_file = os.path.join(self.output_dir, ABUNDANT_COUNTS_FILE)
        # Fragment size distribution
        self.library_metrics_file = os.path.join(self.output_dir, LIBRARY_METRICS_FILE)
        self.frag_size_dist_plot_file = os.path.join(self.output_dir, FRAG_SIZE_DIST_PLOT_FILE)
        # repeat element results
        self.repeat_elements_file = os.path.join(self.output_dir, REPEAT_ELEMENTS_FILE)
        # tophat results
        self.tophat_dir = os.path.join(self.output_dir, TOPHAT_DIR)
        self.tophat_bam_file = os.path.join(self.output_dir, TOPHAT_BAM_FILE)
        self.tophat_bam_index_file = os.path.join(self.output_dir, TOPHAT_BAM_INDEX_FILE)
        self.tophat_juncs_file = os.path.join(self.output_dir, TOPHAT_JUNCTIONS_FILE)
        # Picard metrics
        self.alignment_summary_metrics = os.path.join(self.output_dir, PICARD_ALIGNMENT_SUMMARY_METRICS)
        self.insert_size_histogram_pdf = os.path.join(self.output_dir, PICARD_INSERT_SIZE_HISTOGRAM_PDF)
        self.insert_size_metrics = os.path.join(self.output_dir, PICARD_INSERT_SIZE_METRICS)
        self.quality_by_cycle_metrics = os.path.join(self.output_dir, PICARD_QUALITY_BY_CYCLE_METRICS)
        self.quality_by_cycle_pdf = os.path.join(self.output_dir, PICARD_QUALITY_BY_CYCLE_PDF)
        self.quality_distribution_metrics = os.path.join(self.output_dir, PICARD_QUALITY_DISTRIBUTION_METRICS)
        self.quality_distribution_pdf = os.path.join(self.output_dir, PICARD_QUALITY_DISTRIBUTION_PDF)
        self.rnaseq_metrics = os.path.join(self.output_dir, PICARD_RNASEQ_METRICS)
        self.rnaseq_metrics_pdf = os.path.join(self.output_dir, PICARD_RNASEQ_METRICS_PLOT_PDF)
        # bigwig file prefix
        self.coverage_bigwig_prefix = os.path.join(self.output_dir, COVERAGE_BIGWIG_PREFIX)
        self.coverage_track_file = os.path.join(self.output_dir, COVERAGE_UCSC_TRACK_FILE)        
        # tophat fusion results
        self.tophat_fusion_dir = os.path.join(self.output_dir, TOPHAT_FUSION_DIR)
        self.tophat_fusion_bam_file = os.path.join(self.output_dir, TOPHAT_FUSION_BAM_FILE)        
        self.tophat_fusion_bam_index_file = os.path.join(self.output_dir, TOPHAT_FUSION_BAM_INDEX_FILE)
        self.tophat_fusion_file = os.path.join(self.output_dir, TOPHAT_FUSION_FILE)
        self.tophat_fusion_post_result_file = os.path.join(self.output_dir, TOPHAT_FUSION_POST_RESULT_FILE)
        self.tophat_fusion_tmp_files = (os.path.join(self.output_dir, f) for f in TOPHAT_FUSION_TMP_FILES)
        # cufflinks ab initio output files
        self.cufflinks_ab_initio_dir = os.path.join(self.output_dir, CUFFLINKS_AB_INITIO_DIR)
        self.cufflinks_ab_initio_gtf_file = os.path.join(self.cufflinks_ab_initio_dir, CUFFLINKS_TRANSCRIPTS_GTF_FILE)
        # cufflinks known output files
        self.cufflinks_known_dir = os.path.join(self.output_dir, CUFFLINKS_KNOWN_DIR)
        self.cufflinks_known_gtf_file = os.path.join(self.cufflinks_known_dir, CUFFLINKS_TRANSCRIPTS_GTF_FILE)
        self.cufflinks_known_genes_fpkm_file = os.path.join(self.cufflinks_known_dir, CUFFLINKS_GENES_FILE)
        self.cufflinks_known_isoforms_fpkm_file = os.path.join(self.cufflinks_known_dir, CUFFLINKS_ISOFORMS_FILE)
        # TODO: maintain backwards compatibility
        self.cufflinks_gtf_file = os.path.join(self.output_dir, "cufflinks", "transcripts.gtf")
        # htseq-count output files
        self.htseq_count_known_file = os.path.join(self.output_dir, HTSEQ_COUNT_KNOWN_OUTPUT_FILE)
        # variant calling output
        self.tophat_rmdup_bam_file = os.path.join(self.tmp_dir, TOPHAT_RMDUP_BAM_FILE)
        self.duplicate_metrics = os.path.join(self.output_dir, PICARD_DUPLICATE_METRICS)
        self.varscan_snv_file = os.path.join(self.output_dir, VARSCAN_SNV_FILE)
        # job finished file
        self.job_done_file = os.path.join(self.output_dir, JOB_DONE_FILE)

    def validate(self):
        is_valid = True
        missing_files = []
        # check job done file
        if not os.path.exists(self.job_done_file):
            logging.error("Library %s missing job done file" % (self.library_id))
            missing_files.append(self.job_done_file)
            is_valid = False       
        # read config xml file
        if not os.path.exists(self.config_xml_file):
            logging.error("Library %s missing config xml file" % (self.library_id))
            missing_files.append(self.config_xml_file)
            return False, missing_files
        config = PipelineConfig.from_xml(self.config_xml_file)
        # check FASTQC results
        for f in self.fastqc_data_files:
            if not file_exists_and_nz_size(f):
                logging.error("Library %s missing fastqc data file %s" % (self.library_id, f))
                missing_files.append(f)
                is_valid = False
        for f in self.fastqc_report_files:
            if not file_exists_and_nz_size(f):
                logging.error("Library %s missing fastqc report file %s" % (self.library_id, f))
                missing_files.append(f)
                is_valid = False
        # check sorted abundant reads bam file
        if not check_sam_file(self.sorted_abundant_bam_file, isbam=True):
            logging.error("Library %s missing/corrupt abundant reads BAM file" % (self.library_id))
            missing_files.append(self.sorted_abundant_bam_file)
            is_valid = False
        # check abundant counts file
        if not file_exists_and_nz_size(self.abundant_counts_file):
            logging.error("Library %s missing/corrupt abundant counts file" % (self.library_id))
            missing_files.append(self.abundant_counts_file)
            is_valid = False            
        # check fragment size distribution
        has_library_metrics = check_library_metrics(self.library_metrics_file)
        if not has_library_metrics:
            logging.error("Library %s missing/corrupt inspection results" % (self.library_id))
            missing_files.append(self.library_metrics_file)
            is_valid = False
        # check tophat junctions file
        if not file_exists_and_nz_size(self.tophat_juncs_file):
            logging.error("Library %s missing/empty tophat junctions file" % (self.library_id))
            missing_files.append(self.tophat_juncs_file)
            is_valid = False
        # check tophat bam file
        if not check_sam_file(self.tophat_bam_file, isbam=True):
            logging.error("Library %s missing/corrupt tophat BAM file" % (self.library_id))
            missing_files.append(self.tophat_bam_file)
            is_valid = False
        # check bam index
        if not file_exists_and_nz_size(self.tophat_bam_index_file):
            logging.error("Library %s missing tophat BAM index" % (self.library_id))
            missing_files.append(self.tophat_bam_index_file)
            is_valid = False
        # check picard summary metrics
        if not file_exists_and_nz_size(self.alignment_summary_metrics):
            logging.error("Library %s missing picard alignment summary metrics" % (self.library_id))
            missing_files.append(self.alignment_summary_metrics)
            is_valid = False
        if ((len(self.copied_fastq_files) > 2) and 
            (not file_exists_and_nz_size(self.insert_size_metrics))):
            logging.error("Library %s missing picard insert size metrics" % (self.library_id))
            missing_files.append(self.insert_size_metrics)
            is_valid = False
        if not file_exists_and_nz_size(self.quality_by_cycle_metrics):
            logging.error("Library %s missing picard quality by cycle metrics" % (self.library_id))
            missing_files.append(self.quality_by_cycle_metrics)
            is_valid = False
        if not file_exists_and_nz_size(self.quality_distribution_metrics):
            logging.error("Library %s missing picard quality distribution metrics" % (self.library_id))
            missing_files.append(self.quality_distribution_metrics)
            is_valid = False
        if not file_exists_and_nz_size(self.rnaseq_metrics):
            logging.error("Library %s missing picard rnaseq metrics" % (self.library_id))
            missing_files.append(self.rnaseq_metrics)
            is_valid = False
        # check coverage files
        if not has_library_metrics:
            logging.error("Library %s missing coverage bigwig file(s)" % (self.library_id))
            missing_files.append(self.coverage_bigwig_prefix + ".bw")
            is_valid = False
        else:            
            obj = RnaseqLibraryMetrics.from_file(open(self.library_metrics_file))
            predicted_library_type = obj.predict_library_type(STRAND_SPECIFIC_CUTOFF_FRAC)
            if predicted_library_type == FR_UNSTRANDED:
                bigwig_files = [self.coverage_bigwig_prefix + ".bw"]
            else:
                bigwig_files = [self.coverage_bigwig_prefix + "_pos.bw",
                                self.coverage_bigwig_prefix + "_neg.bw"]
            for f in bigwig_files:
                if not file_exists_and_nz_size(f):
                    logging.error("Library %s missing coverage bigwig file %s" % (self.library_id))
                    is_valid = False
                    missing_files.append(f)
        if not file_exists_and_nz_size(self.coverage_track_file):
            logging.error("Library %s missing coverage track file" % (self.library_id))
            missing_files.append(self.coverage_track_file)
            is_valid = False
        # check tophat fusion (optional)
        if config.tophat_fusion_run:
            # check tophat fusion bam file
            if not check_sam_file(self.tophat_fusion_bam_file, isbam=True):
                logging.error("Library %s missing/corrupt tophat fusion BAM file" % (self.library_id))
                missing_files.append(self.tophat_fusion_bam_file)
                is_valid = False
            # check tophat fusion bam file
            if not file_exists_and_nz_size(self.tophat_fusion_file):
                logging.error("Library %s missing/corrupt tophat fusion file" % (self.library_id))
                missing_files.append(self.tophat_fusion_file)
                is_valid = False
            # check tophat fusion post result file
            if config.tophat_fusion_post_run:
                if not file_exists_and_nz_size(self.tophat_fusion_post_result_file):
                    logging.error("Library %s missing/corrupt tophat fusion post result file" % (self.library_id))
                    missing_files.append(self.tophat_fusion_post_result_file)
                    is_valid = False
        # check cufflinks ab initio files
        if config.cufflinks_ab_initio_run:
            if not file_exists_and_nz_size(self.cufflinks_ab_initio_gtf_file):
                logging.error("Library %s missing cufflinks ab initio gtf file" % (self.library_id))
                missing_files.append(self.cufflinks_ab_initio_gtf_file)
                is_valid = False
        # check cufflinks known files
        if config.cufflinks_known_run:
            if not file_exists_and_nz_size(self.cufflinks_known_gtf_file):
                logging.error("Library %s missing cufflinks known gtf file" % (self.library_id))
                missing_files.append(self.cufflinks_known_gtf_file)
                is_valid = False
        # check htseq files
        if config.htseq_count_run:
            if not file_exists_and_nz_size(self.htseq_count_known_file):
                logging.error("Library %s missing htseq-count output file" % (self.library_id))
                missing_files.append(self.htseq_count_known_file)
                is_valid = False
        # check variant calling files
        if config.varscan_run:
            if not file_exists_and_nz_size(self.duplicate_metrics):
                logging.error("Library %s missing picard duplicate metrics file" % (self.library_id))
                missing_files.append(self.duplicate_metrics)
                is_valid = False
            # check varscan files
            if not file_exists_and_nz_size(self.varscan_snv_file):
                logging.error("Library %s missing varscan snv file" % (self.library_id))
                missing_files.append(self.varscan_snv_file)
                is_valid = False
        return is_valid, missing_files

class GenomeConfig(object):
    fields = ("abundant_bowtie2_index",
              "genome_fasta_file",
              "genome_lexicographical_fasta_file",
              "genome_bowtie1_index",
              "genome_bowtie2_index",
              "fragment_size_bowtie1_index",
              "repbase_bowtie2_index",
              "gene_annotation_refflat",
              "gene_annotation_refgene",
              "gene_annotation_ensgene",
              "picard_ribosomal_intervals",
              "chrom_sizes",
              "known_genes_gtf",
              "transcriptome_bowtie1_index",
              "transcriptome_bowtie2_index",
              "cufflinks_mask_genes")
    
    @staticmethod
    def from_xml_elem(elem):
        g = GenomeConfig()
        g.name = elem.get("name")
        g.root_dir = elem.get("root_dir")
        for attrname in GenomeConfig.fields:
            setattr(g, attrname, elem.findtext(attrname))
        return g

    def to_xml(self, root):
        root.set("name", self.name)
        root.set("root_dir", self.root_dir)
        for attrname in GenomeConfig.fields:
            elem = etree.SubElement(root, attrname)
            elem.text = str(getattr(self, attrname))            

    def resolve_paths(self, root_dir=''):
        # copy genome config object
        g = GenomeConfig()
        g.name = self.name
        g.root_dir = self.root_dir
        # expand paths
        for attrname in GenomeConfig.fields:
            abspath = os.path.join(str(root_dir), g.root_dir, getattr(self, attrname))
            setattr(g, attrname, abspath)
        return g

    def resolve_arg(self, arg):
        newtxt = str(arg)
        for k in re.findall(r'\$\w+', arg):
            v = getattr(self, k[1:])
            newtxt = re.sub('\%s' % k, v, newtxt, count=1)
        return newtxt

    def is_valid(self, references_dir=""):
        valid = True
        abs_root_dir = os.path.join(references_dir, self.root_dir)
        if not os.path.exists(abs_root_dir):
            logging.error("genome root directory %s not found" % (self.root_dir))
            valid = False            
        if not os.path.exists(os.path.join(abs_root_dir, self.abundant_bowtie2_index + ".1.bt2")):
            logging.error("Abundant bowtie2 index %s not found" % (self.abundant_bowtie2_index))
            valid = False
        if not os.path.exists(os.path.join(abs_root_dir, self.abundant_bowtie2_index + ".fa")):
            logging.error("Abundant sequences fasta file %s not found" % (self.abundant_bowtie2_index))
            valid = False
        if not os.path.exists(os.path.join(abs_root_dir, self.genome_bowtie1_index + ".1.ebwt")):
            logging.error("Genome bowtie index %s not found" % (self.genome_bowtie1_index))
            valid = False
        if not os.path.exists(os.path.join(abs_root_dir, self.genome_bowtie2_index + ".1.bt2")):
            logging.error("Genome bowtie2 index %s not found" % (self.genome_bowtie2_index))
            valid = False
        if not os.path.exists(os.path.join(abs_root_dir, self.fragment_size_bowtie1_index + ".1.ebwt")):
            logging.error("Fragment size bowtie index %s not found" % (self.fragment_size_bowtie1_index))
            valid = False
        if not os.path.exists(os.path.join(abs_root_dir, self.repbase_bowtie2_index + ".1.bt2")):
            logging.error("Repeat element bowtie2 index %s not found" % (self.repbase_bowtie2_index))
            valid = False
        if not os.path.exists(os.path.join(abs_root_dir, self.transcriptome_bowtie1_index + ".1.ebwt")):
            logging.error("Transcriptome bowtie index %s not found" % (self.transcriptome_bowtie1_index))
            valid = False
        if not os.path.exists(os.path.join(abs_root_dir, self.transcriptome_bowtie2_index + ".1.bt2")):
            logging.error("Transcriptome bowtie2 index %s not found" % (self.transcriptome_bowtie2_index))
            valid = False
        for attrname in ("genome_fasta_file",
                         "genome_lexicographical_fasta_file",
                         "gene_annotation_refflat",
                         "gene_annotation_refgene",
                         "gene_annotation_ensgene",
                         "picard_ribosomal_intervals",
                         "chrom_sizes",
                         "known_genes_gtf",
                         "cufflinks_mask_genes"):
            if not os.path.exists(os.path.join(abs_root_dir, getattr(self, attrname))):
                logging.error("Annotation file '%s' not found" % (getattr(self, attrname)))
                valid = False
        return valid

class ServerConfig(object):
    @staticmethod
    def from_xml_elem(elem):
        c = ServerConfig()
        c.name = elem.get("name")
        c.modules_init_script = elem.findtext("modules_init_script")
        c.output_dir = elem.findtext("output_dir")
        c.references_dir = elem.findtext("references_dir")
        # seq repository
        c.seq_dirs = collections.defaultdict(lambda: [])
        seq_repos_elem = elem.find("seq_repos")
        if seq_repos_elem is None:
            raise PipelineConfigError("server %s missing 'seq_repo' element" % (c.name))
        for dir_elem in seq_repos_elem.findall("dir"):
            name = dir_elem.get("name")
            c.seq_dirs[name].append(dir_elem.text)
        # pbs support
        c.pbs = False
        c.node_mem = None
        c.node_processors = None
        c.pbs_script_lines = []
        pbs_elem = elem.find("pbs")
        has_pbs = (pbs_elem.get("use") == "yes")
        if has_pbs:
            c.pbs = True
            c.node_mem = float(pbs_elem.findtext("node_mem"))
            c.node_processors = int(pbs_elem.findtext("node_processors"))
            for line_elem in pbs_elem.findall("script_line"):
                c.pbs_script_lines.append(line_elem.text)
        return c

    def to_xml(self, root):
        root.set("name", self.name)
        for attrname in ("modules_init_script",
                         "output_dir",
                         "references_dir"):
            attrval = getattr(self,attrname)
            if attrval is not None:
                elem = etree.SubElement(root, attrname)
                elem.text = str(attrval)
        # seq repository
        seq_repos_elem = etree.SubElement(root, "seq_repos")
        for name, dir_txt in self.seq_dirs.iteritems():
            elem = etree.SubElement(seq_repos_elem, "dir")
            elem.set("name", name)
            elem.text = dir_txt
        # pbs
        pbs_elem = etree.SubElement(root, "pbs")
        use_pbs = "yes" if self.pbs else "no"
        pbs_elem.set("use", use_pbs)
        if self.pbs:
            elem = etree.SubElement(pbs_elem, "node_mem")
            elem.text = self.node_mem
            elem = etree.SubElement(pbs_elem, "node_processors")
            elem.text = self.node_processors            
            for line in self.pbs_script_lines:
                elem = etree.SubElement(pbs_elem, "script_line")
                elem.text = line

    def is_valid(self):
        valid = True
        # check directories
        if os.path.isdir(self.output_dir):
            logging.debug("Checking for 'output' directory %s... found" % (self.output_dir))
        else:
            logging.error("'output' directory %s not found" % (self.output_dir))
            valid = False            
        if os.path.isdir(self.references_dir):
            logging.debug("Checking for 'references' directory %s... found" % (self.references_dir))
        else:
            logging.error("'references' directory %s not found" % (self.references_dir))
            valid = False
        if os.path.exists(self.modules_init_script):
            logging.debug("Checking for modules init script %s... found" % (self.modules_init_script))
        else:
            logging.error("modules init script %s not found" % (self.modules_init_script))
            valid = False
        return valid    

class PipelineConfig(object):
    @staticmethod
    def from_xml(xmlfile):
        tree = etree.parse(xmlfile)  
        root = tree.getroot()
        c = PipelineConfig()
        # genome config
        c.genomes = {}
        for elem in root.findall("genome"):
            g = GenomeConfig.from_xml_elem(elem)
            c.genomes[g.name] = g
        # modules
        modules_elem = root.find("modules")
        c.modules = []
        for elem in modules_elem.findall("module"):
            c.modules.append(elem.text)
        # fragment size parameters
        c.fragment_size_mean_default = int(root.findtext("fragment_size_mean_default"))
        c.fragment_size_stdev_default = int(root.findtext("fragment_size_stdev_default"))
        c.min_fragment_size = int(root.findtext("min_fragment_size"))
        c.max_fragment_size = int(root.findtext("max_fragment_size"))
        # ucsc server url parameters
        c.ucsc_big_data_url = root.findtext("ucsc_big_data_url")
        # tophat parameters
        c.tophat_args = []
        elem = root.find("tophat")
        for arg_elem in elem.findall("arg"):
            c.tophat_args.append(arg_elem.text)
        # tophat fusion parameters
        c.tophat_fusion_args = []
        elem = root.find("tophatfusion")
        c.tophat_fusion_run = parse_bool(elem.get("run", "no"))
        for arg_elem in elem.findall("arg"):
            c.tophat_fusion_args.append(arg_elem.text)
        # tophat-fusion-post parameters
        c.tophat_fusion_post_args = []
        elem = root.find("tophatfusionpost")
        c.tophat_fusion_post_run = parse_bool(elem.get("run", "no"))
        for arg_elem in elem.findall("arg"):
            c.tophat_fusion_post_args.append(arg_elem.text)
        # cufflinks ab initio parameters
        c.cufflinks_ab_initio_args = []
        elem = root.find("cufflinks_ab_initio")
        c.cufflinks_ab_initio_run = parse_bool(elem.get("run", "no"))
        for arg_elem in elem.findall("arg"):
            c.cufflinks_ab_initio_args.append(arg_elem.text)
        # cufflinks known parameters
        c.cufflinks_known_args = []
        elem = root.find("cufflinks_known")
        c.cufflinks_known_run = parse_bool(elem.get("run", "no"))
        for arg_elem in elem.findall("arg"):
            c.cufflinks_known_args.append(arg_elem.text)
        # htseq parameters
        c.htseq_count_args = []
        elem = root.find("htseq")
        c.htseq_count_run = parse_bool(elem.get("run", "no"))
        c.htseq_count_pe = parse_bool(elem.get("pe", "no"))
        for arg_elem in elem.findall("arg"):
            c.htseq_count_args.append(arg_elem.text)
        # varscan parameters
        c.varscan_args = []
        elem = root.find("varscan")
        c.varscan_run = parse_bool(elem.get("run", "no"))
        for arg_elem in elem.findall("arg"):
            c.varscan_args.append(arg_elem.text)
        # server setup
        c.servers = {}
        for elem in root.findall("server"):
            server = ServerConfig.from_xml_elem(elem)
            c.servers[server.name] = server
        return c

    def to_xml(self, output_file):
        root = etree.Element("pipeline")
        # modules
        modules_elem = etree.SubElement(root, "modules")
        for m in self.modules:
            elem = etree.SubElement(modules_elem, "module")
            elem.text = m
        # single element parameters
        for attrname in ("fragment_size_mean_default",
                         "fragment_size_stdev_default",
                         "adaptor_length_default",
                         "min_fragment_size",
                         "max_fragment_size",
                         "ucsc_big_data_url"):
            elem = etree.SubElement(root, attrname)
            elem.text = str(getattr(self, attrname))
        # tophat parameters
        tophat_elem = etree.SubElement(root, "tophat")
        for arg in self.tophat_args:
            elem = etree.SubElement(tophat_elem, "arg")
            elem.text = arg            
        # tophat fusion parameters
        tophat_fusion_elem = etree.SubElement(root, "tophatfusion")
        tophat_fusion_elem.set("run", self.tophat_fusion_run)
        for arg in self.tophat_fusion_args:
            elem = etree.SubElement(tophat_fusion_elem, "arg")
            elem.text = arg
        # tophat fusion post parameters
        tophat_fusion_post_elem = etree.SubElement(root, "tophatfusionpost")
        tophat_fusion_post_elem.set("run", self.tophat_fusion_post_run)
        for arg in self.tophat_fusion_post_args:
            elem = etree.SubElement(tophat_fusion_post_elem, "arg")
            elem.text = arg
        # cufflinks ab initio parameters
        cufflinks_ab_initio_elem = etree.SubElement(root, "cufflinks_ab_initio")
        cufflinks_ab_initio_elem.set("run", self.cufflinks_ab_initio_run)
        for arg in self.cufflinks_ab_initio_args:
            elem = etree.SubElement(cufflinks_ab_initio_elem, "arg")
            elem.text = arg
        # cufflinks known parameters
        cufflinks_known_elem = etree.SubElement(root, "cufflinks_known")
        cufflinks_known_elem.set("run", self.cufflinks_known_run)
        for arg in self.cufflinks_known_args:
            elem = etree.SubElement(cufflinks_known_elem, "arg")
            elem.text = arg
        # htseq parameters
        htseq_elem = etree.SubElement(root, "htseq")
        htseq_elem.set("run", self.htseq_count_run)
        htseq_elem.set("pe", self.htseq_count_pe)
        for arg in self.htseq_count_args:
            elem = etree.SubElement(htseq_elem, "arg")
            elem.text = arg
        # varscan parameters
        varscan_elem = etree.SubElement(root, "varscan")
        varscan_elem.set("run", self.varscan_run)
        for arg in self.varscan_args:
            elem = etree.SubElement(varscan_elem, "arg")
            elem.text = arg
        # servers
        for server in self.servers.itervalues():            
            elem = etree.SubElement(root, "server")
            server.to_xml(elem)
        # genomes
        for genome in self.genomes.itervalues():
            elem = etree.SubElement(root, "genome")
            genome.to_xml(elem)
        # output files
        f = open(output_file, "w")
        # indent for pretty printing
        indent_xml(root)
        print >>f, etree.tostring(root)
        f.close() 

    def is_valid(self, server_name):
        """
        ensure configuration is valid
        """
        valid = True
        # check server
        if server_name not in self.servers:
            logging.error("Server %s not found" % (server_name))
            return False
        server = self.servers[server_name]
        if not server.is_valid():
            logging.error("Server %s missing required paths" % (server_name))
            valid = False
        # check genomes
        for name,genome in self.genomes.iteritems():
            if not genome.is_valid(server.references_dir):
                logging.error("Genome %s missing required files" % (name))
                valid = False
        #
        # Check software installation
        #
        # check fastqc
        msg = 'fastqc'
        if check_executable('fastqc'):
            logging.debug("Checking for '%s' binary... found" % msg)
        else:
            logging.error("'%s' binary not found or not executable" % msg)
            valid = False
        # check samtools
        msg = 'samtools'
        if check_executable('samtools'):
            logging.debug("Checking for '%s' binary... found" % msg)
        else:
            logging.error("'%s' binary not found or not executable" % msg)
            valid = False
        # java
        msg = 'java'
        if check_executable("java"):
            logging.debug("Checking for '%s' binary... found" % msg)
        else:
            logging.error("'%s' binary not found or not executable" % msg)
            valid = False
        # picard
        if "PICARDPATH" not in os.environ:
            logging.debug("PICARDPATH environment variable not set")
            valid = False
        else:
            picard_dir = os.environ["PICARDPATH"]
            jarfile = os.path.join(picard_dir, "CollectMultipleMetrics.jar")
            if os.path.exists(jarfile):
                logging.debug("Checking for picard tools... found")
            else:
                logging.error("Picard jarfile '%s' not found" % (jarfile))
                valid = False
        # varscan
        if "VARSCANPATH" not in os.environ:
            logging.debug("VARSCANPATH environment variable not set")
            valid = False
        else:
            varscan_dir = os.environ["VARSCANPATH"]
            jarfile = os.path.join(varscan_dir, "VarScan.jar")
            if os.path.exists(jarfile):
                logging.debug("Checking for varscan... found")
            else:
                logging.error("VarScan jarfile '%s' not found" % (jarfile))
                valid = False
        # check R
        msg = 'R'
        if check_executable('R'):
            logging.debug("Checking for '%s' binary... found" % msg)
        else:
            logging.error("'%s' binary not found or not executable" % msg)
            valid = False
        msg = 'Rscript'
        if check_executable("Rscript"):
            logging.debug("Checking for '%s' binary... found" % msg)
        else:
            logging.error("'%s' binary not found or not executable" % msg)
            valid = False                
        msg = 'bowtie'
        if check_executable('bowtie'):
            logging.debug("Checking for '%s' binary... found" % msg)
        else:
            logging.error("'%s' binary not found or not executable" % msg)
            valid = False
        msg = 'bowtie2'
        if check_executable('bowtie2'):
            logging.debug("Checking for '%s' binary... found" % msg)
        else:
            logging.error("'%s' binary not found or not executable" % msg)
            valid = False
        msg = 'tophat'
        if check_executable('tophat'):
            logging.debug("Checking for '%s' binary... found" % msg)
        else:
            logging.error("'%s' binary not found or not executable" % msg)
            valid = False
        msg = 'cufflinks'
        if check_executable('cufflinks'):
            logging.debug("Checking for '%s' binary... found" % msg)
        else:
            logging.error("'%s' binary not found or not executable" % msg)
            valid = False
        # check BEDTools
        msg = 'BEDTools'
        if check_executable("bedtools"):
            logging.debug("Checking for '%s' binary... found" % msg)
        else:
            logging.error("'%s' binary not found or not executable" % msg)
            valid = False            
        # check UCSC binary
        msg = 'UCSC bedGraphToBigWig'
        if check_executable('bedGraphToBigWig'):
            logging.debug("Checking for '%s' binaries... found" % msg)
        else:
            logging.error("'%s' binaries not found or not executable" % msg)
            valid = False
        # check htseq-count
        msg = 'htseq-count'
        if check_executable("htseq-count"):
            logging.debug("Checking for '%s' binary... found" % msg)
        else:
            logging.error("'%s' binary not found or not executable" % msg)
            valid = False
        # check for bx python library
        try:
            import bx.intervals.intersection
            logging.debug("Checking for 'bx python' library... found")
        except ImportError, e:
            logging.error("Package 'bx python' not found")
            valid = False
        # check for pysam library
        try:
            import pysam
            logging.debug("Checking for 'pysam' library... found")
        except ImportError, e:
            logging.error("Package 'pysam' not found")
            valid = False
        return valid
