'''
Created on Dec 19, 2012

@author: mkiyer
'''
import sys
import os
import logging
import argparse

from assemblyline.lib.librarytable import LibraryInfo
import assemblyline.rnaseq.lib.picard as picard

def sum_htseq_count(filename):
    count = 0
    for line in open(filename):
        fields = line.strip().split('\t')
        count += int(fields[1])
    return count

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("library_table")
    parser.add_argument("runs_dir")
    parser.add_argument("htseq_dir")
    args = parser.parse_args()
    for lib in LibraryInfo.from_file(args.library_table):
        if not lib.is_valid():
            logging.warning("\tLibrary %s not valid" % (lib.library_id)) 
            continue
        # total picard
        metrics_file = os.path.join(args.runs_dir, lib.library_id, "picard.alignment_summary_metrics")
        if not os.path.exists(metrics_file):
            logging.error("\tLibrary %s no metrics file '%s'" % (lib.library_id, metrics_file))
            continue
        obj = picard.AlignmentSummaryMetrics(metrics_file)
        picard_count = obj.get_total_reads()
        # total htseq
        output_file = os.path.join(args.htseq_dir, lib.library_id, "htseq.txt")
        if not os.path.exists(output_file):
            logging.error("\tLibrary %s no htseq file '%s'" % (lib.library_id, output_file))
            continue
        htseq_count = sum_htseq_count(output_file)
        print '\t'.join(map(str, [lib.library_id, picard_count, htseq_count]))

if __name__ == '__main__':
    sys.exit(main())