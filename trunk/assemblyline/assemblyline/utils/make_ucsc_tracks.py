'''
Created on Feb 6, 2013

@author: mkiyer
'''
import subprocess
import argparse
import sys
import os
import collections

from assemblyline.lib.base import check_executable

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--baseurl", dest="baseurl", default=None) 
    parser.add_argument("prefix")
    parser.add_argument("chrom_sizes_file")
    args = parser.parse_args()
    if not check_executable("bedToBigBed"):
        parser.error("bedToBigBed binary not found in PATH")
    if not check_executable("bedGraphToBigWig"):
        parser.error("'bedGraphToBigWig' executable not found in PATH")
    # find input files
    if not os.path.exists(args.chrom_sizes_file):
        parser.error("chrom sizes file %s not found" % (args.chrom_sizes_file))
    bed_file = args.prefix + ".bed"
    if not os.path.exists(bed_file):
        parser.error("BED file %s not found" % (bed_file))
    bedgraph_files = [args.prefix + "_none.bedgraph",
                      args.prefix + "_neg.bedgraph",
                      args.prefix + "_pos.bedgraph"]
    for bedgraph_file in bedgraph_files:
        if not os.path.exists(bedgraph_file):
            parser.error("Bedgraph file %s not found" % (bedgraph_file))
    # convert to bigbed
    bigbed_file = args.prefix + ".bb"
    retcode = subprocess.call(["bedToBigBed", bed_file, args.chrom_sizes_file, bigbed_file])
    if retcode != 0:
        print >>sys.stderr, "bedToBigBed ERROR"
        return 1
    # print track lines
    track_file = bed_file + ".ucsc_track"
    f = open(track_file)
    fields = f.next().strip().split()
    f.close()
    track_options = ["track"]
    has_type = False
    for i,field in enumerate(fields):
        if field.startswith("type"):
            fields[i] = "type=bigBed"
            has_type = True
            break
    if not has_type:
        track_options.append("type=bigBed")
    track_options.extend(fields[1:])
    if args.baseurl is not None:
        track_options.append('bigDataUrl="%s%s"' % (args.baseurl, os.path.abspath(bigbed_file)))
    track_line = ' '.join(track_options)
    print track_line        
    # convert to bigwig
    for bedgraph_file in bedgraph_files:
        bwfile = os.path.splitext(bedgraph_file)[0] + ".bw"
        retcode = subprocess.call(["bedGraphToBigWig", bedgraph_file, args.chrom_sizes_file, bwfile])
        if retcode != 0:
            print >>sys.stderr, "bedGraphToBigWig ERROR"
            return 1
        track_file = bedgraph_file + ".ucsc_track"
        f = open(track_file)
        fields = f.next().strip().split()
        f.close()
        track_options = ["track"]
        for field in fields[1:]:
            k,v = field.split("=")
            if k == "type":
                track_options.append("type=bigWig")
            else:
                track_options.append(field)
        if args.baseurl is not None:
            track_options.append('bigDataUrl="%s%s"' % (args.baseurl, os.path.abspath(bwfile)))
        track_line = ' '.join(track_options)
        print track_line        
    return 0

if __name__ == '__main__':
    sys.exit(main())
