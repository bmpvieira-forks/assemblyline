'''
Created on Nov 2, 2011

@author: mkiyer
'''

SANGER_ENCODING = "Sanger / Illumina 1.9"
SOLEXA_ENCODING = "Illumina < 1.3"
ILLUMINA_13_ENCODING = "Illumina 1.3"
ILLUMINA_15_ENCODING = "Illumina 1.5"

ENCODING_VALUES = (SANGER_ENCODING, 
                   SOLEXA_ENCODING, 
                   ILLUMINA_13_ENCODING, 
                   ILLUMINA_15_ENCODING)
ENCODING_TO_QUAL_FORMAT = {SANGER_ENCODING: "sanger",
                           SOLEXA_ENCODING: "solexa",
                           ILLUMINA_13_ENCODING: "illumina",
                           ILLUMINA_15_ENCODING: "illumina"}

def get_read_length(fastqc_data_file):
    for line in open(fastqc_data_file):
        if not line: continue
        line = line.strip()
        if line.startswith("Sequence length"):
            return int(line.split()[-1])

def get_total_reads(fastqc_data_file):
    for line in open(fastqc_data_file):
        if not line: continue
        line = line.strip()
        if line.startswith("Total Sequences"):
            return int(line.split()[-1])

def get_fastq_encoding(fastqc_data_file):
    for line in open(fastqc_data_file):
        if not line: continue
        line = line.strip()
        if line.startswith("Encoding"):
            return line.split("\t")[-1]
