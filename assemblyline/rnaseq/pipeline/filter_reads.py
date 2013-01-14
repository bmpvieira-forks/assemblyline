'''
Created on Jan 10, 2013

@author: mkiyer
'''
import argparse
import logging
import os
import sys
import subprocess

# project imports
from assemblyline.rnaseq.lib.base import BOWTIE2_QUAL_MAP, \
    file_exists_and_nz_size, get_fastqc_encoding, \
    ENCODING_TO_QUAL_FORMAT, check_executable
import assemblyline.rnaseq.lib.config as config
import assemblyline.rnaseq.pipeline
_pipeline_dir = assemblyline.rnaseq.pipeline.__path__[0]

def get_quality_score_encoding(fastqc_data_files):
    # get fastq format from fastqc output
    encodings = set()
    for f in fastqc_data_files:
        if file_exists_and_nz_size(f):
            encoding = get_fastqc_encoding(f)
            if encoding not in ENCODING_TO_QUAL_FORMAT:
                logging.error("Unrecognized FASTQ encoding %s" % (encoding))
                return config.JOB_ERROR, None
            encodings.add(ENCODING_TO_QUAL_FORMAT[encoding])
    if len(encodings) > 1:
        logging.error("Detected different FASTQ encodings in paired-end files")
        return config.JOB_ERROR, None
    elif len(encodings) == 0:
        logging.error("Could not locate FASTQC data to determine encoding")
        return config.JOB_ERROR, None
    return config.JOB_SUCCESS, encodings.pop()

def filter_reads(bowtie2_index,
                 fastqc_data_files,
                 input_files,
                 output_files,
                 bam_file,
                 counts_file,
                 picard_dir,
                 tmp_dir,
                 num_processors):
    # get quality scores
    retcode, quals = get_quality_score_encoding(fastqc_data_files)
    if retcode != 0:
        logging.error("Error reading quality score format from FASTQC output")
        return config.JOB_ERROR
    bowtie2_qual_arg = BOWTIE2_QUAL_MAP[quals]
    # map reads using bowtie2
    bowtie_threads = max(1, num_processors - 1)
    logging.info("Mapping reads using bowtie2")
    args = ["bowtie2", "-p", str(bowtie_threads), 
            bowtie2_qual_arg, "-I", "0", "-X", "1000",
            "--end-to-end", "--sensitive", "--reorder",
            "-x", bowtie2_index]
    if len(input_files) == 1:
        args.extend(['-U', input_files[0]])
    else:
        args.extend(['-1', input_files[0],
                     '-2', input_files[1]])
    aln_p = subprocess.Popen(args, stdout=subprocess.PIPE)
    # filter reads
    logging.info("Filtering reads")
    prefix = os.path.basename(os.path.splitext(bam_file)[0])
    unsorted_bam_file = os.path.join(tmp_dir, prefix + ".unsorted.bam") 
    args = ["python", os.path.join(_pipeline_dir, "filter_reads_child.py"),
            "-", unsorted_bam_file, counts_file, ','.join(output_files)]
    retcode1 = subprocess.call(args, stdin=aln_p.stdout)
    retcode2 = aln_p.wait()
    if retcode1 + retcode2 > 0:
        logging.error("filter_reads_child.py failed with error code %d" % (retcode))
        if os.path.exists(unsorted_bam_file):
            os.remove(unsorted_bam_file)
        return 1
    # sort abundant reads bam file
    msg = "Sorting BAM file"
    logging.info(msg)
    args = ["java", "-Xmx4g",
            "-jar", os.path.join(picard_dir, "SortSam.jar"),
            "INPUT=%s" % (unsorted_bam_file),
            "OUTPUT=%s" % (bam_file),
            "SO=coordinate",            
            "CREATE_INDEX=true",
            "MAX_RECORDS_IN_RAM=2000000",
            "TMP_DIR=%s" % tmp_dir]
    command = ' '.join(map(str, args))
    logging.debug("Command string: %s" % (command))
    retcode = subprocess.call(command, shell=True)
    if retcode != 0:
        logging.error("%s FAILED" % (msg))
        if os.path.exists(bam_file):
            os.remove(bam_file)
        return 1
    # remove temp files
    if os.path.exists(unsorted_bam_file):
        os.remove(unsorted_bam_file)
    return 0

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("--num-processors", dest="num_processors", 
                        type=int, default=1)
    parser.add_argument("--picard-dir", dest="picard_dir",
                        default=None)
    parser.add_argument("--tmp-dir", dest="tmp_dir", default=None)
    parser.add_argument("bowtie2_index")
    parser.add_argument("fastqc_data_files")
    parser.add_argument("input_files")
    parser.add_argument("output_files")
    parser.add_argument("bam_file")
    parser.add_argument("counts_file")
    args = parser.parse_args()
    # check command line args
    if not check_executable("bowtie2"):
        parser.error("bowtie2 binary not found")
    if args.picard_dir is None:
        if "PICARDPATH" in os.environ:
            picard_dir = os.environ["PICARDPATH"]
        else:
            parser.error("PICARDPATH not found in environment")
    else:
        if not os.path.exists(args.picard_dir):
            parser.error("Picard directory %s not found" % (args.picard_dir))
        if not os.path.exists(os.path.join(args.picard_dir, "SortSam.jar")):
            parser.error("Picard SortSam.jar file not found")
        picard_dir = args.picard_dir
    if args.tmp_dir is None:
        tmp_dir = os.getcwd()
    else:
        tmp_dir = args.tmp_dir
    input_files = args.input_files.split(',')
    for f in input_files:
        if not os.path.exists(f):
            parser.error("Input file %s not found" % (f))
    fastqc_data_files = args.fastqc_data_files.split(',')
    for f in fastqc_data_files:
        if not os.path.exists(f):
            parser.error("FASTQC data file %s not found" % (f))
    output_files = args.output_files.split(',')
    if len(input_files) > 2:
        parser.error("Can only specify up to 2 input files")
    if len(input_files) != len(fastqc_data_files):
        parser.error("Number of input files != number of FASTQC data files")
    if len(input_files) != len(output_files):
        parser.error("Number of input files != number of output data files")
    # start filtering
    return filter_reads(args.bowtie2_index,
                        fastqc_data_files,
                        input_files,
                        output_files,
                        args.bam_file,
                        args.counts_file,
                        picard_dir,
                        tmp_dir,
                        args.num_processors)

if __name__ == '__main__':
    sys.exit(main())