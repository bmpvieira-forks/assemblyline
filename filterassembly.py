'''
Created on Dec 17, 2010

@author: mkiyer
'''
import collections
import itertools
import argparse
import tempfile
import os
import logging

from lib.bx.intersection import Interval, IntervalTree
import lib.bed as bed

def build_interval_trees(bedfile):
    trees = collections.defaultdict(lambda: IntervalTree())
    exons = collections.defaultdict(lambda: set())
    for g in bed.BEDFeature.parse(open(bedfile)):
        for start, end in g.exons:
            if (start, end) in exons[g.chrom]:
                continue
            exons[g.chrom].add((start, end))    
            trees[g.chrom].insert_interval(Interval(start, end))
    return trees

#def get_overlap_fraction(a, b):
#    '''
#    a,b are BEDFeature objects
#    '''
#    # insert b into an interval tree
#    introntree = IntervalTree()
#    for x,y in itertools.izip(b.exons[:-1], b.exons[1:]): 
#        introntree.insert_interval(Interval(x[1], y[0], strand=b.strand))
#    # a and b must be intron compatible
#    if any((len(introntree.find(start, end)) > 0) for start, end in a.exons):
#        return 0.0
#     check exons
#    btree = IntervalTree()
#    for bstart, bend in b.exons:
#        btree.insert_interval(Interval(bstart, bend, strand=b.strand))
#     search for matches to each exon
#    overlap_bases = 0
#    tx_length = 0    
#    for astart, aend in a.exons:
#        bintervals = btree.find(astart, aend)
#        max_exon_overlap_bases = 0
#        for interval in bintervals:            
#            bstart, bend = interval.start, interval.end
#            overlap_start = astart if bstart < astart else bstart
#            overlap_end = bend if bend < aend else aend
#            max_exon_overlap_bases = max(max_exon_overlap_bases, (overlap_end - overlap_start))
#        overlap_bases += max_exon_overlap_bases
#        tx_length += (aend - astart)
#    return overlap_bases / float(tx_length)

def get_overlap_fraction(g, trees):
    # Determine percent overlap with annotated genes
    overlap_fraction = 0.0
    if g.chrom in trees:
        tree = trees[g.chrom]            
        overlap_bases = 0
        tx_length = 0    
        for astart, aend in g.exons:                
            bintervals = tree.find(astart, aend)
            max_exon_overlap_bases = 0
            for interval in bintervals:            
                bstart, bend = interval.start, interval.end
                overlap_start = astart if bstart < astart else bstart
                overlap_end = bend if bend < aend else aend
                max_exon_overlap_bases = max(max_exon_overlap_bases, (overlap_end - overlap_start))
            overlap_bases += max_exon_overlap_bases
            tx_length += (aend - astart)            
        overlap_fraction = overlap_bases / float(tx_length)
    return overlap_fraction

def filter_assembly(bedfile, annotation_bed_file,
                    overlap_fraction_threshold):    
    annotation_trees = build_interval_trees(annotation_bed_file)    
    # make temp file to store transcript data
    #fd,filename = tempfile.mkstemp(suffix='.txt', prefix='transcript_stats',
    #                               dir=os.path.dirname(bedfile))
    #os.close(fd)
    #fileh = open(filename, 'w')
    print '\t'.join(['name', 'length', 'num_exons', 'score', 'lanes', 'ann', 'overlap'])
    # write transcript data
    for g in bed.BEDFeature.parse(open(bedfile)):
        # determine percent overlap with annotated genes
        overlap_fraction = get_overlap_fraction(g, annotation_trees)
        annotated = '1' if (overlap_fraction >= overlap_fraction_threshold) else '0'
        # get transcript length
        length = sum((end - start) for start, end in g.exons)
        # print transcript data to temporary file
        print '\t'.join([g.name, 
                         str(length),
                         str(g.exon_count), 
                         str(g.score),
                         str(g.attr_fields[1]),                         
                         annotated, 
                         str(round(overlap_fraction, 3))])

#R_SCRIPT = os.path.join(os.path.dirname(__file__), "classify_juncs.R")

#def classify_juncs(infile, outfile, epsfile, summaryfile):
#    fd,tmpfile = tempfile.mkstemp(suffix=".pred", prefix="tmp", dir=os.path.dirname(outfile))
#    os.close(fd)
#    logging.debug("Running classification tree on %s" % (infile))
#    args = ['R', '--slave', '--vanilla', '--args', infile, tmpfile, epsfile]
#    retcode = subprocess.call(args, stdin=open(R_SCRIPT), stdout=open(summaryfile, 'w'))        

        #print annotated, overlap_fraction, g
    #logger.debug("%s: Running classification tree" % (task_name))
    #args = ['R', '--slave', '--vanilla', '--args', datafile, resultfile]
    #retcode = subprocess.call(args, stdin=open(r_script), stdout=open(summaryfile, 'w'))
    

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")    
    parser = argparse.ArgumentParser()
    parser.add_argument("--overlap-threshold", dest="overlap_fraction_threshold", type=float, default=0.5)
    parser.add_argument("assembly_bed_file")
    parser.add_argument("annotation_bed_file")
    options = parser.parse_args()
    filter_assembly(options.assembly_bed_file, 
                    options.annotation_bed_file,
                    options.overlap_fraction_threshold)

if __name__ == '__main__': main()
