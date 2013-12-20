'''
Created on Sep 11, 2013

@author: mkiyer
'''
import re
import os
import sys
import argparse
import logging
import glob
import subprocess
import itertools
import numpy as np
from multiprocessing import Process, Queue

BIGWIG_TO_BEDGRAPH_BIN = 'bigWigToBedGraph'

def which(program):
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)
    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file
    return None

class BEDFeature(object):    
    __slots__ = ('chrom', 'tx_start', 'tx_end', 'name', 'score', 'strand',
                 'cds_start', 'cds_end', 'exon_count', 'block_starts', 
                 'block_sizes', 'exons', 'attr_fields')

    def __str__(self):
        fields = [self.chrom,
                  str(self.tx_start),
                  str(self.tx_end),
                  self.name,
                  str(self.score),
                  self.strand,
                  str(self.cds_start),
                  str(self.cds_end),
                  '0',
                  str(self.exon_count),
                  ','.join(map(str, self.block_sizes)) + ',',
                  ','.join(map(str, self.block_starts)) + ',']
        return '\t'.join(fields)

    @staticmethod
    def from_string(line):        
        if line is None:
            return None
        line = line.strip()
        if line.startswith('#'):
            logging.debug("skipping comment line: %s" % (line))
            return None
        if line.startswith('track'):
            logging.debug("skipping track header line: %s"  % (line))
            return None
        fields = line.split('\t')
        # first six fields are required
        g = BEDFeature()
        g.chrom = fields[0]
        g.tx_start = int(fields[1])
        g.tx_end = int(fields[2])
        if len(fields) <= 3:
            g.name = '%s:%d-%d' % (g.chrom, g.tx_start, g.tx_end)
        else:
            g.name = fields[3]
        if len(fields) <= 4:
            g.score = 0
            g.strand = '.'
        else:
            g.score = fields[4]
            g.strand = fields[5]        
        if len(fields) <= 6:
            g.cds_start = g.tx_start
            g.cds_end = g.tx_end
            g.exon_count = 1
            g.exons = [(g.tx_start, g.tx_end)]
        else:
            g.cds_start = int(fields[6])
            g.cds_end = int(fields[7])
            g.exon_count = int(fields[9])
            g.block_sizes = map(int, fields[10].split(',')[:-1])
            g.block_starts = map(int, fields[11].split(',')[:-1])            
            g.exons = []
            for start, size in itertools.izip(g.block_starts, g.block_sizes):
                g.exons.append((g.tx_start + start, g.tx_start + start + size))
        if len(fields) <= 12:
            g.attr_fields = []
        else:
            g.attr_fields = fields[12:]
        return g

    @staticmethod
    def parse(line_iter):
        for line in line_iter:
            if not line:
                continue
            if not line.strip():
                continue            
            if line.startswith("#"):
                continue
            if line.startswith("track"):
                continue
            yield BEDFeature.from_string(line)

def find_bigwig_files(root_path, pattern):            
    # find bigwig files
    glob_pattern = pattern.replace("${CHROM}", "*")
    re_pattern = pattern.replace("${CHROM}", r'(.+)')
    chrom_bigwig_dict = {}
    for filename in glob.glob(os.path.join(root_path, glob_pattern)):
        m = re.match(re_pattern, filename)
        if m is not None:
            chrom = os.path.basename(m.group(1))
            chrom_bigwig_dict[chrom] = filename
    return chrom_bigwig_dict

def best_sliding_window(arr, window_size, score_func='np.mean'):
    if window_size > len(arr):
        return 'NA'
    if len(arr) == 0:
        return 'NA'
    scores = []
    num_windows = len(arr) - window_size + 1 
    for i in xrange(num_windows):
        scores.append(score_func(arr[i:i+window_size]))
    return np.max(scores)

def extract_bigwig_data(feature, bigwig_file):
    exon_arrays = []
    for exon_start,exon_end in feature.exons:
        exon_size = exon_end - exon_start
        arr = np.zeros(exon_size, dtype=np.float)
        args = [BIGWIG_TO_BEDGRAPH_BIN, 
                '-chrom=%s' % (feature.chrom),
                '-start=%d' % (exon_start),
                '-end=%d' % (exon_end),
                bigwig_file,
                '/dev/stdout']
        p = subprocess.Popen(args, stdout=subprocess.PIPE)
        res = p.communicate()[0]
        for line in res.split('\n'):
            if not line:
                continue
            fields = line.strip().split('\t')
            line_start = int(fields[1])
            line_end = int(fields[2])
            value = float(fields[3])
            arr_start = line_start - exon_start
            arr_end = line_end - exon_start
            arr[arr_start:arr_end] = value
        if feature.strand == '-':
            arr = arr[::-1]
        exon_arrays.append(arr)
    if feature.strand == '-':
        exon_arrays.reverse()
    # compute summary statistics
    arr = np.concatenate(exon_arrays)
    return arr

def conservation_parallel(bed_file, window_sizes, chrom_bigwig_dict, 
                          bigwig_to_bedgraph_bin, num_processes):
    def _producer(q):
        for line in open(bed_file):
            q.put(line.strip())
        for i in xrange(num_processes):
            q.put(None)
    def _consumer(input_queue, output_queue):
        while True:
            line = input_queue.get()
            if line is None:
                break
            f = BEDFeature.from_string(line)
            # retrieve conservation data
            bigwig_file = chrom_bigwig_dict[f.chrom]
            arr = extract_bigwig_data(f, bigwig_file, bigwig_to_bedgraph_bin)
            # measure conservation at various sliding windows
            window_scores = []
            for window_size in window_sizes:
                window_scores.append(best_sliding_window(arr, window_size, np.mean))
            # measure average conservation
            fields = [f.name, '%s:%d-%d[%s]' % (f.chrom, f.tx_start, f.tx_end, f.strand),
                      str(len(arr)), str(np.mean(arr))]
            fields.extend(map(str,window_scores))
            result = '\t'.join(fields)
            output_queue.put(result)
        output_queue.put(None)
    # create multiprocessing queues for passing data
    input_queue = Queue(maxsize=num_processes*3)
    output_queue = Queue(maxsize=num_processes*3) 
    # start producer process
    logging.debug("Starting producer process and %d workers" % (num_processes))
    producer = Process(target=_producer, args=(input_queue,))
    producer.start()    
    # start consumer processes
    procs = []
    for i in xrange(num_processes):
        p = Process(target=_consumer, args=(input_queue, output_queue))
        p.start()
        procs.append(p)
    # output header
    fields = ['name', 'position', 'transcript_length', 'mean']
    fields.extend(map(str,window_sizes))
    print '\t'.join(fields)
    # get results from consumers
    num_alive = num_processes
    while num_alive > 0:
        result = output_queue.get()
        if result is None:
            num_alive -= 1
            logging.debug("Main process detected worker finished, %d still alive" % (num_alive))
        else:
            print result
    logging.debug("Joining all processes")
    # wait for producer to finish
    producer.join()
    # wait for consumers to finish
    for p in procs:
        p.join()

def conservation_serial(bed_file, window_sizes, chrom_bigwig_dict):
    # output header fields
    fields = ['name', 'position', 'transcript_length', 'mean']
    fields.extend(map(str,window_sizes))
    print '\t'.join(fields)
    # process bed file
    for f in BEDFeature.parse(open(bed_file)):
        # retrieve conservation data
        bigwig_file = chrom_bigwig_dict[f.chrom]
        arr = extract_bigwig_data(f, bigwig_file)
        # measure conservation at various sliding windows
        window_scores = []
        for window_size in window_sizes:
            window_scores.append(best_sliding_window(arr, window_size, np.mean))
        # measure average conservation
        fields = [f.name, '%s:%d-%d[%s]' % (f.chrom, f.tx_start, f.tx_end, f.strand),
                  str(len(arr)), str(np.mean(arr))]
        fields.extend(map(str,window_scores))
        print '\t'.join(fields)    

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", dest="num_processes", type=int, default=1)
    parser.add_argument("--pattern", dest="pattern", default=r'${CHROM}.phastCons46way.bw')
    parser.add_argument('--window-sizes', dest='window_sizes', default='30,60,90,150,300,600,900,1500,3000')
    parser.add_argument('--bed', dest='bed_file', default=None)
    parser.add_argument("bigwig_file_dir")
    args = parser.parse_args()
    if which(BIGWIG_TO_BEDGRAPH_BIN) is None:
        parser.error('bigWigToBedGraph binary not found in PATH')
    if args.bed_file is not None:
        if not os.path.exists(args.bed_file):
            parser.error('bed file %s not found' % (args.bed_file))
    else:
        parser.error('specify a bed file using --bed')
    window_sizes = map(int, args.window_sizes.split(','))
    num_processes = max(1, args.num_processes)
    # find bigwig files
    chrom_bigwig_dict = find_bigwig_files(args.bigwig_file_dir, args.pattern)
    if num_processes > 1:
        conservation_parallel(args.bed_file, window_sizes, chrom_bigwig_dict, num_processes)
    else:
        conservation_serial(args.bed_file, window_sizes, chrom_bigwig_dict)
    return 0

if __name__ == '__main__':
    sys.exit(main())