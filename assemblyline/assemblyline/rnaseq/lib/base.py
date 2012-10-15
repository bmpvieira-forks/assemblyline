'''
Created on Aug 3, 2011

@author: mkiyer
'''
import subprocess
import os
import gzip
import bz2
import zipfile
from math import log10
from string import maketrans

# Quality score formats
SANGER_FORMAT = "sanger"
SOLEXA_FORMAT = "solexa"
ILLUMINA_FORMAT = "illumina"
FASTQ_QUAL_FORMATS = [SANGER_FORMAT, SOLEXA_FORMAT, ILLUMINA_FORMAT]
BOWTIE2_QUAL_MAP = {SANGER_FORMAT: '--phred33',
                    SOLEXA_FORMAT: '--solexa-quals',
                    ILLUMINA_FORMAT: '--phred64'}

# FASTQC constants
SANGER_ENCODING = "Sanger / Illumina 1.9"
SOLEXA_ENCODING = "Illumina < 1.3"
ILLUMINA_13_ENCODING = "Illumina 1.3"
ILLUMINA_15_ENCODING = "Illumina 1.5"
ENCODING_VALUES = (SANGER_ENCODING, 
                   SOLEXA_ENCODING, 
                   ILLUMINA_13_ENCODING, 
                   ILLUMINA_15_ENCODING)
ENCODING_TO_QUAL_FORMAT = {SANGER_ENCODING: SANGER_FORMAT,
                           SOLEXA_ENCODING: SOLEXA_FORMAT,
                           ILLUMINA_13_ENCODING: ILLUMINA_FORMAT,
                           ILLUMINA_15_ENCODING: ILLUMINA_FORMAT}

#Translation table for reverse Complement, with ambiguity codes
DNA_COMPLEMENT = maketrans( "ACGTRYKMBDHVacgtrykmbdhv", "TGCAYRMKVHDBtgcayrmkvhdb" )
RNA_COMPLEMENT = maketrans( "ACGURYKMBDHVacgurykmbdhv", "UGCAYRMKVHDBugcayrmkvhdb" )
#Translation table for DNA <--> RNA
DNA_TO_RNA = maketrans( "Tt", "Uu" )
RNA_TO_DNA = maketrans( "Uu", "Tt" )

# FASTQ chars per line
BASES_PER_LINE = 50

def DNA_complement( sequence ):
    '''complement DNA sequence string'''
    return sequence.translate( DNA_COMPLEMENT )
def DNA_reverse_complement( sequence ):
    '''returns the reverse complement of the sequence'''
    return DNA_complement(sequence[::-1])
def to_DNA( sequence ):
    return sequence.translate( DNA_TO_RNA )
#complement RNA sequence string
def RNA_complement( sequence ):
    return sequence.translate( RNA_COMPLEMENT )
def RNA_reverse_complement( self, sequence ):
    return RNA_complement( sequence[::-1] )
def to_RNA( sequence ):
    return sequence.translate( RNA_TO_DNA )

def split_seq(seq, chars_per_line=BASES_PER_LINE):
    pos = 0
    newseq = []
    while pos < len(seq):
        if pos + chars_per_line > len(seq):        
            endpos = len(seq)
        else:
            endpos = pos + chars_per_line
        newseq.append(seq[pos:endpos])
        pos = endpos
    return '\n'.join(newseq)

def get_solexa_qual_conversion_table():
    """
    return a translation table that can be used by str.translate() for
    converting solexa to sanger quality scores
    """
    offset = 64
    conv_table = ['!'] * 256
    conv_table[offset:] = "I" * (256-offset)
    for solq in xrange(-5, 40):
        phredq = 10*log10(1 + 10**(solq/10.0))
        phredchr = chr(int(round(33 + phredq)))
        conv_table[offset + solq] = phredchr
    conv_string = ''.join(conv_table)
    return maketrans(''.join(map(chr, range(256))), conv_string)

def get_illumina_qual_conversion_table():
    """Illumina 1.3+ format can encode a Phred quality score from 0 to 62 
    using ASCII 64 to 126 (although in raw read data Phred scores from 0 
    to 40 only are expected).
    """
    offset = 64
    conv_table = ['!'] * 256
    for x in xrange(0, 62):
        conv_table[offset+x] = chr(33 + x)
    conv_table[offset+40:] = "I" * (256-(offset+40))
    conv_string = ''.join(conv_table)
    return maketrans(''.join(map(chr, range(256))), conv_string)    

def get_sanger_qual_conversion_table():
    offset = 33
    tbl = map(chr, range(256))
    tbl[:offset] = "!" * offset
    tbl[offset+40:] = "I" * (256-(offset+40))
    return maketrans(''.join(map(chr, range(256))), ''.join(tbl))

def get_qual_conversion_func(qual_format):
    conv_tables = {SANGER_FORMAT: get_sanger_qual_conversion_table(),
                   ILLUMINA_FORMAT: get_illumina_qual_conversion_table(),
                   SOLEXA_FORMAT: get_solexa_qual_conversion_table()}
    tbl = conv_tables[qual_format]
    return lambda q: q.translate(tbl)

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

def detect_format(f):
    if f.endswith(".gz"):
        return "gz"
    elif f.endswith(".bz2"):
        return "bz2"
    elif f.endswith(".zip"):
        return "zip"
    else:
        return "txt"

def open_compressed(fastq_file):
    compression_format = detect_format(fastq_file)
    if compression_format == "gz":
        filehandle = gzip.open(fastq_file, "r")
    elif compression_format == "bz2":
        filehandle = bz2.BZ2File(fastq_file, "r")
    elif compression_format == "zip":
        filehandle = zipfile.ZipFile(fastq_file, "r")
    else:
        filehandle = open(fastq_file, "r")
    return filehandle

def detect_read_length(fastq_file):
    f = open_compressed(fastq_file)
    f.next()
    seq = f.next()
    f.close()
    return len(seq)

def check_sam_file(filename, isbam=False):
    is_valid = True
    if not file_exists_and_nz_size(filename):
        is_valid = False
    else:
        import pysam
        try:
            fmt = "rb" if isbam else "r"
            samfh = pysam.Samfile(filename, fmt)
            samfh.close()   
        except:
            is_valid = False
    return is_valid

def parse_reads_by_qname(samfh):
    """
    generator function to parse and return lists of
    reads that share the same qname
    """    
    reads = []
    for read in samfh:        
        if len(reads) > 0 and read.qname != reads[-1].qname:
            yield reads
            reads = []
        reads.append(read)
    if len(reads) > 0:
        yield reads

def parse_pe_reads(bamfh):
    pe_reads = ([], [])
    # reads must be sorted by qname
    num_reads = 0
    prev_qname = None
    for read in bamfh:
        # get read attributes
        qname = read.qname
        readnum = 1 if read.is_read2 else 0
        # if query name changes we have completely finished
        # the fragment and can reset the read data
        if num_reads > 0 and qname != prev_qname:
            yield pe_reads
            # reset state variables
            pe_reads = ([], [])
            num_reads = 0
        pe_reads[readnum].append(read)
        prev_qname = qname
        num_reads += 1
    if num_reads > 0:
        yield pe_reads

def check_executable(filename):
    # check that samtools binary exists
    devnullfh = open(os.devnull, 'w')        
    try:
        p = subprocess.Popen([filename], stdout=devnullfh, stderr=devnullfh)
        p.kill()
    except OSError:
        return False
    devnullfh.close()
    return True

def up_to_date(outfile, infile, nzsize=True):
    if not os.path.exists(infile):
        return False
    if not os.path.exists(outfile):
        return False
    if nzsize and (os.path.getsize(outfile) == 0):
        return False
    return os.path.getmtime(outfile) >= os.path.getmtime(infile)

def many_up_to_date(output_files, input_files):
    skip = True
    for input_file in input_files:
        for output_file in output_files:
            skip = skip and up_to_date(output_file, input_file)
    return skip

def file_exists_and_nz_size(filename):
    if filename is None:
        return False
    if not os.path.exists(filename):
        return False
    if os.path.getsize(filename) == 0:
        return False
    return True

def parse_bool(s):    
    return True if s[0].lower() == "t" else False

def parse_string_none(s):
    return None if s == "None" else s

# in-place XML prettyprint formatter
def indent_xml(elem, level=0):
    i = "\n" + level*"  "
    if len(elem):
        if not elem.text or not elem.text.strip():
            elem.text = i + "  "
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
        for elem in elem:
            indent_xml(elem, level+1)
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
    else:
        if level and (not elem.tail or not elem.tail.strip()):
            elem.tail = i