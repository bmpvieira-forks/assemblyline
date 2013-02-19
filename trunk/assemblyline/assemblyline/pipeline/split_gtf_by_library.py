'''
Created on Nov 30, 2011

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
import logging
import argparse
import os
import sys
import collections
import multiprocessing

import assemblyline
import assemblyline.lib.config as config
from assemblyline.lib.gtf import GTFFeature
from assemblyline.lib.base import CategoryCounts, GTFAttr

class FileHandleCache(object):
    def __init__(self, keyfunc, maxsize=config.MAX_OPEN_FILE_DESCRIPTORS):
        self.fileh_dict = collections.OrderedDict()
        self.file_dict = {}
        self.keyfunc = keyfunc
        self.maxsize = maxsize
        self.hits = 0
        self.misses = 0

    def get_file_handle(self, k):
        if k in self.fileh_dict:
            fileh = self.fileh_dict.pop(k)
            self.hits += 1
        else:
            if k not in self.file_dict:
                filename = self.keyfunc(k)
                self.file_dict[k] = filename 
                mode = 'w'
            else:
                # not opening for the first time so append
                filename = self.file_dict[k]
                mode = 'a'
                self.misses += 1
            # control number of open files
            if len(self.fileh_dict) > self.maxsize:
                # close least recently accessed file
                lrufileh = self.fileh_dict.popitem(0)[1]
                lrufileh.close()        
            # open file                
            fileh = open(filename, mode)
        # update
        self.fileh_dict[k] = fileh
        return fileh
    
    def close(self):
        for fileh in self.fileh_dict.itervalues():
            fileh.close()

def split_gtf_file(gtf_file, split_dir, ref_gtf_file, category_counts_file):
    # split input gtf by library and mark test ids
    keyfunc = lambda myid: os.path.join(split_dir, "%s.gtf" % (myid))
    cache = FileHandleCache(keyfunc)
    ref_fileh = open(ref_gtf_file, 'w')
    counts_dict = collections.defaultdict(lambda: CategoryCounts())
    logging.info("Splitting transcripts by library")
    for f in GTFFeature.parse(open(gtf_file)):
        is_ref = bool(int(f.attrs[GTFAttr.REF]))
        if is_ref:
            print >>ref_fileh, str(f)
            continue
        library_id = f.attrs[GTFAttr.LIBRARY_ID]
        # keep statistics
        if f.feature_type == 'transcript':
            category = int(f.attrs[GTFAttr.CATEGORY])           
            countsobj = counts_dict[library_id]
            countsobj.library_id = library_id
            countsobj.counts[category] += 1
        # write features from each library to separate files
        fileh = cache.get_file_handle(library_id)
        print >>fileh, str(f)
    # close open file handles
    ref_fileh.close()
    cache.close()
    logging.debug("File handle cache hits: %d" % (cache.hits))
    logging.debug("File handle cache misses: %d" % (cache.misses))
    # write library category statistics
    logging.info("Writing category count statistics")
    fh = open(category_counts_file, "w")
    print >>fh, '\t'.join(CategoryCounts.header_fields())
    for countsobj in counts_dict.itervalues():
        fields = countsobj.to_fields()
        print >>fh, '\t'.join(map(str, fields))
    fh.close()
       
def main():
    multiprocessing.freeze_support()
    # parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("run_dir")
    parser.add_argument("-v", "--verbose", action="store_true", 
                        dest="verbose", default=False)
    args = parser.parse_args()
    # check command line parameters
    if not os.path.exists(args.run_dir):
        parser.error("Run directory %s not found" % (args.run_dir))
    # set logging level
    if args.verbose:
        level = logging.DEBUG
    else:
        level = logging.INFO
    logging.basicConfig(level=level,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    logging.info("AssemblyLine %s" % (assemblyline.__version__))
    logging.info("----------------------------------")   
    # show parameters
    logging.info("Parameters:")
    logging.info("run directory:    %s" % (args.run_dir))
    logging.info("verbose logging:  %s" % (args.verbose))
    logging.info("----------------------------------")   
    # setup results
    results = config.AssemblylineResults(args.run_dir)
    if not os.path.exists(results.classify_dir):
        os.makedirs(results.classify_dir)
    # split gtf file
    split_gtf_file(results.annotated_transcripts_gtf_file, 
                   results.classify_dir,
                   results.ref_gtf_file,
                   results.category_counts_file)
    logging.info("Done")
    return 0

if __name__ == "__main__":
    sys.exit(main())
