'''
Created on Sep 11, 2013

@author: mkiyer
'''
import os
import sys
import argparse
import logging
import numpy as np
from multiprocessing import Process, Queue

from assemblyline.utils.conservation.base import BIGWIG_TO_BEDGRAPH_BIN, \
    BEDFeature, which, extract_bigwig_data, find_bigwig_files

def best_sliding_window(arr, window_size, score_func='np.mean'):
    if window_size > len(arr):
        return 'NA'
    if len(arr) == 0:
        return 'NA'
    scores = []
    num_windows = len(arr) - window_size + 1
    for i in xrange(num_windows):
        warr = arr[i:i+window_size]
        if np.any(np.isnan(warr)):
            continue
        scores.append(score_func(arr[i:i+window_size]))
    if len(scores) == 0:
        return np.nan
    return np.max(scores)

def conservation_parallel(bed_file, window_sizes, chrom_bigwig_dict, 
                          num_processes):
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
            arr = extract_bigwig_data(f, bigwig_file)
            # measure conservation at various sliding windows
            window_scores = []
            for window_size in window_sizes:
                window_scores.append(best_sliding_window(arr, window_size, np.mean))
            # measure average conservation
            finitearr = arr[np.isfinite(arr)]
            if len(finitearr) == 0:
                mean_cons = np.nan
            else:
                mean_cons = np.mean(finitearr)
            fields = [f.name, '%s:%d-%d[%s]' % (f.chrom, f.tx_start, f.tx_end, f.strand),
                      str(len(arr)), str(mean_cons)]
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
        # calc mean conservation
        finitearr = arr[np.isfinite(arr)]
        if len(finitearr) == 0:
            mean_cons = np.nan
        else:
            mean_cons = np.mean(finitearr)
        fields = [f.name, '%s:%d-%d[%s]' % (f.chrom, f.tx_start, f.tx_end, f.strand),
                  str(len(arr)), str(mean_cons)]
        fields.extend(map(str,window_scores))
        print '\t'.join(fields)    

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", dest="num_processes", type=int, default=1)
    parser.add_argument("--pattern", dest="pattern", default=r'{{CHROM}}.phastCons46way.bw')
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