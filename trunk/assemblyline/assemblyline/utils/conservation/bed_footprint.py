'''
Created on Jan 10, 2014

@author: mkiyer
'''
import argparse
import logging
import collections

from assemblyline.lib.bx.cluster import ClusterTree
from base import BEDFeature

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument('bed_file')
    args = parser.parse_args()

    trees = collections.defaultdict(lambda: ClusterTree(0,1))    
    for f in BEDFeature.parse(open(args.bed_file)):
        tree = trees[f.chrom]        
        for start,end in f.exons:
            tree.insert(start, end, 1)
    footprint = 0
    for chrom in sorted(trees):
        chromprint = 0
        tree = trees[chrom]
        for start,end,indexes in tree.getregions():
            chromprint += (end - start)
        #print chrom, chromprint
        footprint += chromprint
    print 'total', footprint
        
                
if __name__ == '__main__':
    main()     