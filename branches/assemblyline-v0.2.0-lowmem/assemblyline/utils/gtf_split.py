'''
Created on Sep 6, 2012

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
import argparse
import logging
import os
import sys

from assemblyline.lib.gtf import GTFFeature

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("--gtf-split-attr", dest="gtf_split_attr", default="library_id")
    parser.add_argument("gtf_file")
    args = parser.parse_args()
    fhdict = {}
    gtf_split_attr = args.gtf_split_attr
    for f in GTFFeature.parse(open(args.gtf_file)):
        if gtf_split_attr not in f.attrs:
            val = "na_missing"
        else:
            val = f.attrs[args.gtf_split_attr]
        filename = "%s.gtf" % (val)
        if not os.path.exists(filename):
            fh = open(filename, "w")
            fhdict[val] = fh
        else:
            fh = fhdict[val]
        print >>fh, str(f)
    for fh in fhdict.itervalues():
        fh.close()
    return 0

if __name__ == '__main__':
    sys.exit(main())