'''
Created on Dec 4, 2012

@author: mkiyer
'''
import logging
import os
import argparse

LOG_EVERY = 1000000

def parse_lines(line_iter, numlines=1):
    """
    generator that returns list of 'numlines' lines at a time
    """
    try:
        while True:
            yield [line_iter.next().rstrip() for x in xrange(numlines)]
    except StopIteration:
        pass

def check_fastq(fastq_files):
    # setup file iterators for input fastq files
    infhs = [open(f) for f in fastq_files]
    fqiters = [parse_lines(f, numlines=4) for f in infhs]
    linenum = 0
    log_next = LOG_EVERY
    try:
        while True:
            pelines = [it.next() for it in fqiters]
            # get read1 first line of fq record, and remove "@" symbol
            read1_name = pelines[0][0][1:]
            read2_name = pelines[1][0][1:]
            # remove whitespace and/or read number tags /1 or /2
            read1_name = read1_name.split()[0].split("/")[0]
            read2_name = read2_name.split()[0].split("/")[0]
            if read1_name != read2_name:
                logging.error("FASTQ entry %d (line %d) has reads that do not match: %s != %s" % 
                              (linenum, linenum*4, read1_name, read2_name))
            linenum += 1
            if linenum == log_next:
                logging.debug("Read %d records" % (linenum))
                log_next += LOG_EVERY
    except StopIteration:
        pass
    # cleanup
    for fh in infhs:
        fh.close()
    logging.debug("Inspected %d fragments" % (linenum))
    return 0

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")    
    parser = argparse.ArgumentParser()
    parser.add_argument("fastq_files", nargs="+")    
    args = parser.parse_args()
    return check_fastq(args.fastq_files)

if __name__ == '__main__':
    main()