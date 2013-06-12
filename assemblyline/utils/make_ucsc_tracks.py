'''
Created on Feb 6, 2013

@author: mkiyer
'''
import subprocess
import argparse
import sys
import os

from assemblyline.lib.base import check_executable

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--baseurl", dest="baseurl", default=None) 
    parser.add_argument("output_dir")
    parser.add_argument("chrom_sizes_file")
    args = parser.parse_args()
    if not check_executable("bedToBigBed"):
        parser.error("bedToBigBed binary not found in PATH")
    if not check_executable("bedGraphToBigWig"):
        parser.error("'bedGraphToBigWig' executable not found in PATH")
    if not os.path.exists(args.chrom_sizes_file):
        parser.error("chrom sizes file %s not found" % (args.chrom_sizes_file))
    if not os.path.exists(args.output_dir):
        parser.error("output dir %s not found" % (args.output_dir))
    # input files
    prefix = os.path.dirname(args.output_dir)
    output_dir = os.path.abspath(args.output_dir)
    bed_file = os.path.join(output_dir, "assembly.bed")
    bed_track_file = bed_file + ".ucsc_track"
    bedgraph_files = [os.path.join(output_dir, "assembly_none.bedgraph"),
                      os.path.join(output_dir, "assembly_neg.bedgraph"),
                      os.path.join(output_dir, "assembly_pos.bedgraph")]
    bedgraph_track_files = [x + ".ucsc_track" for x in bedgraph_files]
    if not os.path.exists(bed_file):
        parser.error("BED file %s not found" % (bed_file))
    if not os.path.exists(bed_track_file):
        parser.error("BED track file %s not found" % (bed_track_file))
    for i in xrange(len(bedgraph_files)):
        if not os.path.exists(bedgraph_files[i]):
            parser.error("Bedgraph file %s not found" % (bedgraph_files[i]))
        if not os.path.exists(bedgraph_track_files[i]):
            parser.error("Bedgraph track file %s not found" % (bedgraph_track_files[i]))
    # convert to bigbed
    bigbed_file = os.path.join(output_dir, "assembly.bb")
    retcode = subprocess.call(["bedToBigBed", bed_file, args.chrom_sizes_file, bigbed_file])
    if retcode != 0:
        print >>sys.stderr, "bedToBigBed ERROR"
        return 1
    # print track lines
    f = open(bed_track_file)
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
