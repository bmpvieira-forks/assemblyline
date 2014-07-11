'''
Created on Jan 7, 2014

@author: mkiyer
'''
import argparse
import logging
import os
import sys
import collections
from multiprocessing import Process, Queue
import numpy as np

from assemblyline.utils.conservation.base import BIGWIG_TO_BEDGRAPH_BIN, \
    BEDFeature, which, extract_bigwig_data, find_bigwig_files

# histogram bins
NUM_BINS = 20001
BINS = np.linspace(-10.0, 10.0, NUM_BINS)
BIN_MIN = BINS[0]
BIN_MAX = BINS[-1]

def bed_feature_conservation(f, chrom_bigwig_dict, hists):
    # retrieve conservation data
    bigwig_file = chrom_bigwig_dict[f.chrom]
    arr = extract_bigwig_data(f, bigwig_file)
    # ignore missing values
    finitearr = arr[np.isfinite(arr)]
    if len(finitearr) == 0:
        cons_str = 'NA'
    else:
        hists[f.name] += np.histogram(np.clip(finitearr, BIN_MIN, BIN_MAX), BINS)[0]
        cons_str = str(np.sum(finitearr))
    fields = f.name.split('|')
    fields.extend([f.chrom, str(f.tx_start), str(f.tx_end), f.strand, str(len(finitearr)), cons_str])
    return fields

def conservation_parallel(bed_file, chrom_bigwig_dict, num_processes,
                          results_file, hists_file):
    def _producer(q):
        for line in open(bed_file):
            q.put(line.strip())
        for i in xrange(num_processes):
            q.put(None)
    def _consumer(worker_index, input_queue, output_queue):
        hists = collections.defaultdict(lambda: np.zeros(NUM_BINS-1, dtype=np.float))
        while True:
            line = input_queue.get()
            if line is None:
                break
            f = BEDFeature.from_string(line)
            fields = bed_feature_conservation(f, chrom_bigwig_dict, hists)
            result = '\t'.join(fields)
            output_queue.put(result)
        np.savez('w%d.npz' % (worker_index), **hists)    
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
        p = Process(target=_consumer, args=(i, input_queue, output_queue))
        p.start()
        procs.append(p)
    # get results from consumers
    num_alive = num_processes
    outfile = open(results_file, 'w')
    while num_alive > 0:
        result = output_queue.get()
        if result is None:
            num_alive -= 1
            logging.debug("Main process detected worker finished, %d still alive" % (num_alive))
        else:
            print >>outfile, result
    outfile.close()
    logging.debug("Joining all processes")
    # wait for producer to finish
    producer.join()
    # wait for consumers to finish
    for p in procs:
        p.join()
    # merge histograms
    logging.debug("Merging %d worker histograms" % (num_processes))
    hists = collections.defaultdict(lambda: np.zeros(NUM_BINS-1, dtype=np.float))
    for i in xrange(num_processes):
        # aggregate numpy arrays
        hist_file = 'w%d.npz' % (i)
        npzfile = np.load(hist_file)
        for k in npzfile.iterkeys():
            hists[k] += np.array(npzfile[k])
        npzfile.close()
        os.remove(hist_file)
    hists = dict(hists)
    np.savez(hists_file, **hists)

def main():
    # parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("-v", "--verbose", action="store_true", 
                        dest="verbose", default=False)
    parser.add_argument('-p', '--num-processes', type=int, 
                        dest='num_processes', default=1)
    parser.add_argument("--pattern", dest="pattern", 
                        default=r'{{CHROM}}.phyloP46way.bw')
    parser.add_argument("bigwig_file_dir")
    parser.add_argument("bed_file")
    args = parser.parse_args()
    # set logging level
    if args.verbose:
        level = logging.DEBUG
    else:
        level = logging.INFO
    logging.basicConfig(level=level,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    # check command line parameters
    if which(BIGWIG_TO_BEDGRAPH_BIN) is None:
        parser.error('bigWigToBedGraph binary not found in PATH')
    if not os.path.exists(args.bed_file):
        parser.error("BED file %s not found" % (args.bed_file))
    prefix = os.path.splitext(args.bed_file)[0]
    results_file = prefix + '.results.txt'
    hists_file = prefix + '.hists.npz'
    # find bigwig files
    logging.info("Indexing bigWig files")
    chrom_bigwig_dict = find_bigwig_files(args.bigwig_file_dir, args.pattern)
    # process bed file
    logging.info("Measuring conservation")
    if args.num_processes > 1:
        conservation_parallel(args.bed_file, chrom_bigwig_dict, args.num_processes,
                              results_file, hists_file)
    else:       
        hists = collections.defaultdict(lambda: np.zeros(NUM_BINS-1, dtype=np.float))
        with open(results_file, 'w') as outfile:
            for f in BEDFeature.parse(open(args.bed_file)):
                fields = bed_feature_conservation(f, chrom_bigwig_dict, hists)
                print >>outfile, '\t'.join(fields)
        np.savez(hists_file, **hists)    
    return 0

if __name__ == "__main__":
    sys.exit(main())
