'''
Created on Aug 10, 2011

@author: mkiyer
'''
import argparse
import logging
import os
import subprocess

# local imports
import pysam

# project imports
from assemblyline.rnaseq.lib.base import parse_reads_by_qname, \
    DNA_reverse_complement, BOWTIE2_QUAL_MAP, file_exists_and_nz_size, \
    get_fastq_encoding, ENCODING_TO_QUAL_FORMAT
import assemblyline.rnaseq.lib.config as config

def to_fastq(r, f):
    seq = DNA_reverse_complement(r.seq) if r.is_reverse else r.seq
    qual = r.qual[::-1] if r.is_reverse else r.qual 
    print >>f, "@%s\n%s\n+\n%s" % (r.qname, seq, qual)

def get_fastqc_encoding(fastqc_data_files):
    # get fastq format from fastqc output
    encodings = set()
    for f in fastqc_data_files:
        if file_exists_and_nz_size(f):
            encoding = get_fastq_encoding(f)
            if encoding not in ENCODING_TO_QUAL_FORMAT:
                logging.error("Unrecognized FASTQ encoding %s" % (encoding))
                return config.JOB_ERROR
            encodings.add(ENCODING_TO_QUAL_FORMAT[encoding])
    if len(encodings) > 1:
        logging.error("Detected different FASTQ encodings in paired-end files")
        return config.JOB_ERROR, None
    elif len(encodings) == 0:
        logging.error("Could not locate FASTQC data to determine encoding")
        return config.JOB_ERROR, None
    return config.JOB_SUCCESS, encodings.pop()

def filter_reads(input_files, fastq_files, bam_file):
    input_sam_fhs = []
    fastq_fhs = []
    for readnum in xrange(len(input_files)):
        input_sam_fhs.append(pysam.Samfile(input_files[readnum], "r"))
        fastq_fhs.append(open(fastq_files[readnum], "w"))
    bamfh = pysam.Samfile(bam_file, "wb", template=input_sam_fhs[0])
    read_iters = [parse_reads_by_qname(f) for f in input_sam_fhs]
    try:
        while True:            
            pe_reads = [it.next() for it in read_iters]
            # if either read maps to the abundant sequences, discard the pair
            pe_mapped_reads = []
            for readnum, reads in enumerate(pe_reads):
                mapped_reads = []
                for r in reads:
                    if not r.is_unmapped:
                        r.is_paired = True
                        r.mate_is_unmapped = True
                        if readnum == 1:
                            r.is_read2 = True
                        else:
                            r.is_read1 = True
                        mapped_reads.append(r)
                pe_mapped_reads.append(mapped_reads)
            if all((len(mapped_reads) == 0) for mapped_reads in pe_mapped_reads):
                for readnum, reads in enumerate(pe_reads):
                    to_fastq(reads[0], fastq_fhs[readnum])            
            else:
                for readnum,mapped_reads in enumerate(pe_mapped_reads):
                    if len(mapped_reads) == 0:
                        bamfh.write(pe_reads[readnum][0])
                    else:
                        for r in mapped_reads:
                            bamfh.write(r)
    except StopIteration:
        pass
    bamfh.close()
    for fh in fastq_fhs:
        fh.close()
    for fh in input_sam_fhs:
        fh.close()

def filter_abundant_sequences(fastqc_data_files,
                              input_files,
                              output_files,
                              sorted_abundant_bam_file,
                              abundant_index,
                              tmp_dir,
                              picard_dir,
                              num_processors):
    if "PICARDPATH" in os.environ:
        picard_dir = os.environ["PICARDPATH"]
    #
    # get quality scores
    #
    retcode, quals = get_fastqc_encoding(fastqc_data_files)
    if retcode != 0:
        logging.error("Error reading quality score format from FASTQC output")
        return config.JOB_ERROR
    bowtie2_qual_arg = BOWTIE2_QUAL_MAP[quals]
    # map reads against abundant sequences
    abundant_sam_files = []
    msg = "Mapping reads against abundant sequences using bowtie2"
    for readnum in xrange(len(input_files)):
        abundant_sam_file = os.path.join(tmp_dir, config.ABUNDANT_SAM_FILES[readnum])
        abundant_sam_files.append(abundant_sam_file)
        args = ["bowtie2", "-p", num_processors, bowtie2_qual_arg,
                "--end-to-end", "--sensitive", "--reorder", 
                "-x", abundant_index,
                "-U", input_files[readnum],
                "-S", abundant_sam_file]
        command = ' '.join(map(str, args))
        logging.debug("Command string: %s" % (command))
        retcode = subprocess.call(command, shell=True)
        if retcode != 0:
            logging.error("%s FAILED" % (msg))
            return config.JOB_ERROR
    #
    # filter abundant reads
    #
    abundant_bam_file = os.path.join(tmp_dir, config.ABUNDANT_BAM_FILE)
    filter_reads(abundant_sam_files, output_files, abundant_bam_file)
    # remove temp files
    for f in abundant_sam_files:
        os.remove(f)
    #
    # sort abundant reads bam file
    #
    msg = "Sorting abundant hits BAM file"
    logging.info(msg)
    args = ["java", "-Xmx4g",
            "-jar", os.path.join(picard_dir, "SortSam.jar"),
            "INPUT=%s" % (abundant_bam_file),
            "OUTPUT=%s" % (sorted_abundant_bam_file),
            "SO=coordinate",            
            "CREATE_INDEX=true",
            "MAX_RECORDS_IN_RAM=2000000",
            "TMP_DIR=%s" % tmp_dir]
    command = ' '.join(map(str, args))
    logging.debug("Command string: %s" % (command))
    retcode = subprocess.call(command, shell=True)
    if retcode != 0:
        logging.error("%s FAILED" % (msg))
        return config.JOB_ERROR
    # remove temp files
    os.remove(abundant_bam_file)    

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("--num-processors", dest="num_processors", 
                        type=int, default=1)
    parser.add_argument("--picard-dir", dest="picard_dir",
                        default=None)
    parser.add_argument("fastqc_data_files")
    parser.add_argument("input_files")
    parser.add_argument("output_files")
    parser.add_argument("bam_file")
    parser.add_argument("abundant_index")
    parser.add_argument("tmp_dir")
    args = parser.parse_args()
    if args.picard_dir is None:
        if "PICARDPATH" in os.environ:
            picard_dir = os.environ["PICARDPATH"]
        else:
            parser.error("PICARDPATH not found in environment")
    else:
        picard_dir = args.picard_dir
    filter_abundant_sequences(args.fastqc_data_files.split(","),
                              args.input_files.split(","),
                              args.output_files.split(","),
                              args.bam_file,
                              args.abundant_index,
                              args.tmp_dir,
                              picard_dir,
                              args.num_processors)
    return 0

if __name__ == '__main__':
    main()