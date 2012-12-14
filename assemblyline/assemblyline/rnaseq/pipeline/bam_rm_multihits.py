'''
Created on Dec 4, 2012

@author: mkiyer
'''
import logging
import argparse
import os
import sys

import pysam

def parse_pe_readnum_in_qname(bamfh):
    # reads must be sorted by qname
    pe_reads = ([], [])
    num_reads = 0
    prev_qname = None
    for read in bamfh:
        # use "/1" or "/2" suffix to determine read number
        suffix = read.qname[-2:]
        if (suffix == "/1"):
            readnum = 0
            read.is_read1 = True
            read.is_read2 = False
            read.is_paired = True
        elif (suffix == "/2"):
            readnum = 1
            read.is_read1 = False
            read.is_read2 = True
            read.is_paired = True
        else:
            raise Exception("suffix /1 or /2 not found in qname")
        read.qname = read.qname[:-2]
        # if query name changes we have completely finished
        # the fragment and can reset the read data
        if num_reads > 0 and read.qname != prev_qname:
            yield pe_reads
            # reset state variables
            pe_reads = ([], [])
            num_reads = 0
        pe_reads[readnum].append(read)
        prev_qname = read.qname
        num_reads += 1
    if num_reads > 0:
        yield pe_reads

def parse_pe_readnum_in_flags(bamfh, fix_qname=True):
    # reads must be sorted by qname
    pe_reads = ([], [])
    num_reads = 0
    prev_qname = None
    for read in bamfh:
        # get readnum
        readnum = 1 if read.is_read2 else 0
        # optionally remove suffix
        if fix_qname:
            suffix = read.qname[-2:]
            if (suffix == "/1" or suffix == "/2"):
                read.qname = read.qname[:-2]
        # if query name changes we have completely finished
        # the fragment and can reset the read data
        if num_reads > 0 and read.qname != prev_qname:
            yield pe_reads
            # reset state variables
            pe_reads = ([], [])
            num_reads = 0
        pe_reads[readnum].append(read)
        prev_qname = read.qname
        num_reads += 1
    if num_reads > 0:
        yield pe_reads

def bam_pe_rm_multihits(bamfh, outfh, readnum_in_qname):
    num_frags = 0
    orphan_frags = 0
    if readnum_in_qname:
        parse_pe_reads_iter = parse_pe_readnum_in_qname(bamfh)
    else:
        parse_pe_reads_iter = parse_pe_readnum_in_flags(bamfh)        
    for pe_reads in parse_pe_reads_iter:
        n1 = len(pe_reads[0])
        n2 = len(pe_reads[1])
        if (n1 == 0 or n2 == 0):
            orphan_frags += 1
            continue
        r1 = pe_reads[0][0]
        r2 = pe_reads[1][0]
        r1.is_secondary = False
        r2.is_secondary = False
        outfh.write(r1)
        outfh.write(r2)
        num_frags += 1
    outfh.close()
    bamfh.close()
    logging.info("Found %d fragments" % (num_frags))
    logging.info("Skipped %d orphan fragments" % (orphan_frags))

def parse_sr_reads(bamfh, fix_qname=True):
    reads = []
    for read in bamfh:
        if fix_qname:
            suffix = read.qname[-2:]
            if (suffix == "/1" or suffix == "/2"):
                read.qname = read.qname[:-2]        
        if len(reads) > 0 and read.qname != reads[-1].qname:
            yield reads
            reads = []
        reads.append(read)
    if len(reads) > 0:
        yield reads

def bam_sr_rm_multihits(bamfh, outfh):
    num_frags = 0
    for reads in parse_sr_reads(bamfh):        
        r = reads[0]
        r.is_secondary = False
        outfh.write(r)
        num_frags += 1
    outfh.close()
    bamfh.close()
    logging.debug("Found %d fragments" % (num_frags))

def main():
    # setup logging
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    # command line parsing
    parser = argparse.ArgumentParser()
    parser.add_argument('--fragment-layout', dest="fragment_layout", default="paired")
    parser.add_argument('--readnum-in-qname', dest="readnum_in_qname", action="store_true", default=False)
    parser.add_argument('bam_file')
    parser.add_argument('output_file')
    args = parser.parse_args()
    if not os.path.exists(args.bam_file):
        parser.error("Input bam file not found")
    # open files
    bamfh = pysam.Samfile(args.bam_file, "rb")
    sam = (args.output_file == "-" or 
           os.path.splitext(args.output_file)[-1] == ".sam")
    if sam:
        outfh = pysam.Samfile(args.output_file, "wh", template=bamfh)
    else:
        outfh = pysam.Samfile(args.output_file, "wb", template=bamfh)
    try:
        if args.fragment_layout == "paired":
            bam_pe_rm_multihits(bamfh, outfh, args.readnum_in_qname)
        else:
            bam_sr_rm_multihits(bamfh, outfh)
    except Exception as inst:
        logging.error("removing multihits failed with exception %s" % (str(inst)))
        bamfh.close()
        outfh.close()
        if os.path.exists(args.output_file):
            os.remove(args.output_file)
        return 1
    bamfh.close()
    outfh.close()
    return 0

if __name__ == '__main__':
    sys.exit(main())
