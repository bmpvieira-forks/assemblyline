'''
Created on Oct 8, 2012

@author: mkiyer

AssemblyLine: transcriptome meta-assembly from RNA-Seq

Copyright (C) 2012 Matthew Iyer

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
'''
import subprocess
import argparse
import sys
import os

from assemblyline.lib.base import check_executable

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("bedgraph_file")
    parser.add_argument("chrom_sizes_file")
    args = parser.parse_args()
    # check args
    if not os.path.exists(args.bedgraph_file):
        parser.error("bedgraph file %s not found" % (args.bedgraph_file))
    if not os.path.exists(args.chrom_sizes_file):
        parser.error("chrome sizes file %s not found" % (args.chrom_sizes_file))
    # check for executables
    if not check_executable("bedGraphToBigWig"):
        parser.error("'bedGraphToBigWig' executable not found in PATH")
    if not check_executable("sed"):
        parser.error("'sed' executable not found in PATH")
    if not check_executable("sort"):
        parser.error("'sort' executable not found in PATH")
    # convert bedgraph to bigwig
    prefix = os.path.splitext(args.bedgraph_file)[0]
    tmp_noheader_file = prefix + ".noheader.bedgraph"
    bigwig_file = prefix + ".bw"
    # remove header line of file
    outf = open(tmp_noheader_file, "w")
    retcode1 = subprocess.call(["sed", "1,1d", tmp_noheader_file], stdout=outf)
    outf.close()
    # convert to bigwig
    retcode2 = subprocess.call(["bedGraphToBigWig", tmp_noheader_file, args.chrom_sizes_file, bigwig_file])
    # cleanup
    os.remove(tmp_noheader_file)
    if (retcode1 != 0) or (retcode2 != 0):
        return 1
    return 0

if __name__ == '__main__':
    sys.exit(main())
