'''
Created on Aug 4, 2011

@author: mkiyer
'''
from math import log10
from string import maketrans
import gzip
import bz2
import zipfile

#Translation table for reverse Complement, with ambiguity codes
DNA_COMPLEMENT = maketrans( "ACGTRYKMBDHVacgtrykmbdhv", "TGCAYRMKVHDBtgcayrmkvhdb" )
RNA_COMPLEMENT = maketrans( "ACGURYKMBDHVacgurykmbdhv", "UGCAYRMKVHDBugcayrmkvhdb" )
#Translation table for DNA <--> RNA
DNA_TO_RNA = maketrans( "Tt", "Uu" )
RNA_TO_DNA = maketrans( "Uu", "Tt" )

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

BASES_PER_LINE = 50
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

# Quality score formats
SANGER_FORMAT = "sanger"
SOLEXA_FORMAT = "solexa"
ILLUMINA_FORMAT = "illumina"
FASTQ_QUAL_FORMATS = [SANGER_FORMAT, SOLEXA_FORMAT, ILLUMINA_FORMAT]

BOWTIE2_QUAL_MAP = {SANGER_FORMAT: '--phred33',
                    SOLEXA_FORMAT: '--solexa-quals',
                    ILLUMINA_FORMAT: '--phred64'}


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