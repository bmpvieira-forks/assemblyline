'''
Created on Apr 30, 2010

@author: mkiyer
'''
import logging
import os
import sys
import collections
import xml.etree.cElementTree as etree
import multiprocessing
import numpy as np
import shutil

# packages
from bx.intervals.intersection import Interval, IntervalTree
from bx.intervals.cluster import ClusterTree

# project imports
import sequel.genome as genome
from sequel.trackfactory import TrackFactory
#from sequel.io.gtf import write_gtf2_feature
import sequel.io.bed as bedio 
from sequel.io.samplepools import read_sample_pools

# local imports
import pipelineconfig
import cuffcompare
from transcript_parser import parse_gtf_loci

def parse_config_file(xml_file):
    params = {}
    tree = etree.parse(xml_file)
    root = tree.getroot()
    # general parameters
    params['genome_version'] = root.findtext('genome_version')
    logging.info("Genome version: %s" % params['genome_version'])    

    # sample pools
    sample_pools_file = root.findtext('samples')
    sample_pools = read_sample_pools(open(sample_pools_file))
    ordered_samples = []
    logging.info("Ordered samples from cuffcompare:")
    for sample_name, library_ids in zip(sample_pools.names, sample_pools.library_ids):        
        for library_id in library_ids:
            ordered_samples.append(sample_name)
            logging.info("Library %s -> Sample %s" % (library_id, sample_name))
    params['ordered_samples'] = ordered_samples

    # cuffcompare output path
    elem = root.find('cuffcompare')
    params['output_path'] = pipelineconfig.fix_path(elem.findtext('output_path'))    
    logging.info("Cuffcompare output path: %s" % params['output_path'])
    # mapability
    melem = root.find('mapability')
    params['mapability_file'] = melem.findtext('file')
    params['mapability_track'] = melem.findtext('track')
    params['mapability_min'] = float(melem.findtext('min'))
    params['percentile_fpkm'] = float(melem.findtext('percentile_fpkm'))
    logging.info("Mapability file: %s track: %s min mapability: %f percentile rpkm: %f" % 
                 (params['mapability_file'], 
                  params['mapability_track'], 
                  params['mapability_min'],
                  params['percentile_fpkm']))
    data_dir = root.findtext('data_dir')
    annotations = []
    for elem in root.find('annotations'):
        filename = os.path.join(data_dir, elem.text)
        if not os.path.exists(filename):
            logging.error("Could not find annotation file %s" % (filename))
            continue
        annotations.append((elem.tag, filename))
    params['annotations'] = annotations
    logging.info("Annotations: %s" % (annotations))    
    # output directory for results
    output_dir = root.findtext('cuffcompare_analysis_dir')
    if not os.path.exists(output_dir):
        logging.info("Creating cuffcompare analysis dir %s" % output_dir)
        os.makedirs(output_dir)
    params['output_dir'] = output_dir
    return params

def build_interval_trees(trees, bedfile):
    for g in bedio.parse_bed12_file(open(bedfile)):
        for e in g.exons:
            trees[g.chrom].insert_interval(Interval(e[0], e[1]))

def build_annotation_interval_trees(annotations):
    # setup annotation interval trees
    trees = collections.defaultdict(lambda: IntervalTree())
    for name, bedfile in annotations:
        logging.debug("Adding annotation file %s (%s) to interval trees" % (name, bedfile))
        build_interval_trees(trees, bedfile)
    return trees

def gen_tracking_data(cuffcompare_output_path, chrom, ordered_sample_names):
    # read input file to make library_id -> index map
    # get set of recurrent transcripts from tracking file
    tracking_file = os.path.join(cuffcompare_output_path, chrom + '.tracking')
    for tracking_data in cuffcompare.parse_tracking_file(open(tracking_file), ordered_sample_names):
        yield tracking_data

def gen_overlapping_features(features):    
    # break locus into clusters of overlapping features
    cluster_tree = ClusterTree(0,1)     
    for id, f in enumerate(features):
        cluster_tree.insert(f.start, f.end, id)
    # iterate through clusters
    compatstrand = {'+': ['.'],
                    '-': ['.'],
                    '.': ['+', '-']}
    for start, end, ids in cluster_tree.getregions():
        strand_index_map = collections.defaultdict(lambda: [])
        # bin by compatible strands
        for id in ids:
            f = features[id]
            # add this feature id its own strand list
            strand_index_map[f.strand].append(id)
            # add feature to compatible strand lists
            for otherstrand in compatstrand[f.strand]:
                strand_index_map[otherstrand].append(id)
        # yield lists of compatible ids for this cluster
        yield start, end, dict(strand_index_map)
        
def gen_locus_transcripts(cuffcompare_output_path, chrom):
    # break chromosome into loci of overlapping features
    combined_gtf_path = os.path.join(cuffcompare_output_path, chrom + '.combined.gtf')
    for locus_features in parse_gtf_loci(open(combined_gtf_path), 
                                         cuffcompare.cuffcompare_combined_attrs):
        yield locus_features

def measure_uniqueness(chrom, start, end, arrtrack, min_neglogp):
    uarr = arrtrack.fetch(chrom, start, end)
    unique_indexes = np.flatnonzero(uarr >= min_neglogp)        
    longest_stretch = 0
    current_stretch = 0
    current_i = None
    for i in unique_indexes:
        if i-1 == current_i:
            current_stretch += 1
        else:
            if current_stretch > longest_stretch:
                longest_stretch = current_stretch
            current_stretch = 0                    
        current_i = i
    if (current_stretch > longest_stretch):
        longest_stretch = current_stretch
    uavg = np.average(uarr)
    ustdev = np.std(uarr)
    return uavg, ustdev, longest_stretch

def measure_feature(chrom, f,
                    annotation_trees, 
                    tracking_stats, 
                    mapabilitytrack, 
                    mapability_min):   
    # do individual measurements on each feature
    t_id = '%s_%s' % (chrom, f.attrs['transcript_id'])            
    recurrence, fpkm, fpkm_lo, fpkm_hi = tracking_stats.get(t_id, (0, 0, 0, 0))
    uavg, ustdev, ulongest = measure_uniqueness(chrom, f.start, f.end, mapabilitytrack, mapability_min)           
    if chrom in annotation_trees:
        annotated = len(annotation_trees[chrom].find(f.start, f.end)) > 0
    else:
        annotated = False    
    # record feature measurements
    f.attrs['transcript_id'] = t_id
    f.attrs['recur'] = recurrence
    f.attrs['fpkm_lo'] = round(fpkm_lo, 3)
    f.attrs['fpkm_hi'] = round(fpkm_hi, 3)
    f.attrs['fpkm'] = round(fpkm, 3)
    f.attrs['uavg'] = round(uavg, 2)
    f.attrs['ustd'] = round(ustdev,2)
    f.attrs['ulong'] = ulongest
    f.attrs['ann'] = int(annotated)

def measure_cluster_features(strand, feature_indexes, locus_features):
    # measure features of the entire cluster
    max_recurrence = 0
    max_fpkm, max_fpkm_lo, max_fpkm_hi = 0.0, 0.0, 0.0
    max_length = 0
    max_annotated = 0      
    # measure properties using entire list of compatible features
    # which includes the features being studied
    for index in feature_indexes:
        f = locus_features[index]
        # update cluster measurements
        max_recurrence = max(max_recurrence, f.attrs['recur'])
        max_fpkm = max(max_fpkm, f.attrs['fpkm'])
        max_fpkm_lo = max(max_fpkm_lo, f.attrs['fpkm_lo'])
        max_fpkm_hi = max(max_fpkm_hi, f.attrs['fpkm_hi'])
        max_length = max(max_length, f.end - f.start)
        max_annotated = max(max_annotated, f.attrs['ann']) 
    # record cluster measurements for this strand only
    for index in feature_indexes:
        f = locus_features[index]
        if f.strand != strand:
            continue
        # record feature measurements
        f.attrs['max_recur'] = max_recurrence
        f.attrs['max_fpkm_lo'] = round(max_fpkm_lo, 3)
        f.attrs['max_fpkm_hi'] = round(max_fpkm_hi, 3)
        f.attrs['max_fpkm'] = round(max_fpkm, 3)
        f.attrs['max_length'] = max_length
        f.attrs['max_ann'] = max_annotated


def measure_features(chrom,
                     cuffcompare_output_path, 
                     trees, 
                     tracking_stats, 
                     mapabilitytrack, 
                     mapability_min):
    # break into loci
    logging.debug("Reading transcripts on chromosome %s..." % (chrom))    
    for locus_features in gen_locus_transcripts(cuffcompare_output_path, chrom):        
        # features are returned in list of lists (locus -> transcript -> feature)
        # and we need a flat ordered list for clustering
        flat_features = []
        for transcript_features in locus_features:
            flat_features.extend(transcript_features)
            for f in transcript_features:
                # first record all the individual feature measurements
                # measure the feature (measurements are recorded in the feature dictionaries)
                measure_feature(chrom, f, trees, tracking_stats, mapabilitytrack, mapability_min)
        # break locus into clusters of overlapping strand-compatible 
        # transcripts and record measurements of those clusters
        for start, end, strand_index_map in gen_overlapping_features(flat_features):
            for strand, strand_indexes in strand_index_map.iteritems():
                logging.debug("%s:%d-%d(%s) examining %d features" %
                              (chrom, start, end, strand, len(strand_indexes)))
                # measure the features
                measure_cluster_features(strand, strand_indexes, flat_features)
        # return all the features (in the original order)
        yield flat_features


def measure_transcripts(chrom, params, outfileh):
    cuffcompare_output_path = params['output_path']
    mapability_min = params['mapability_min']
    percentile_fpkm = params['percentile_fpkm']
    # setup mapability track
    tf = TrackFactory(params['mapability_file'], 'r')
    mapabilitytrack = tf.get_track(params['mapability_track'])
    # setup annotation interval trees
    trees = build_annotation_interval_trees(params['annotations'])
    # get transcript recurrence from tracking data
    logging.debug("Parsing tracking data for chromosome %s..." % (chrom))    
    tracking_stats = {}
    for tdata in gen_tracking_data(cuffcompare_output_path, chrom, params['ordered_samples']):
        t_id = '%s_%s' % (chrom, tdata.transcript_id)
        ind = np.floor(percentile_fpkm * len(tdata.fpkm))
        fpkm = np.sort(tdata.fpkm)[ind]
        fpkm_lo = np.sort(tdata.fpkm_lo)[ind]
        fpkm_hi = np.sort(tdata.fpkm_hi)[ind]            
        tracking_stats[t_id] = (tdata.recurrence, fpkm, fpkm_lo, fpkm_hi)        
    # read transcripts and measure parameters
    logging.debug("Measuring transcripts on chromosome %s..." % (chrom))    
    for features in measure_features(chrom,
                                     cuffcompare_output_path, 
                                     trees, 
                                     tracking_stats, 
                                     mapabilitytrack, 
                                     mapability_min):
        for f in features:
            print >>outfileh, f
    # cleanup mapability track
    tf.close()

def measure_transcripts_worker(task):
    chrom, params, outfile = task
    outfileh = open(outfile, 'w')
    measure_transcripts(chrom, params, outfileh)
    outfileh.close()
    return outfile

def measure_transcripts_parallel(config_file, outfileh, processors=1, keep_tmp=True):
    # read relevant parts in config file
    params = parse_config_file(config_file)
    # create tasks (one per chromosome)
    genome_version = params['genome_version']
    chrom_names = sorted(genome.get_genome(genome_version).get_chrom_names())    
    #chrom_names = ["chrM"]
    tasklist = []
    for chrom in chrom_names:
        output_file = os.path.join(params["output_dir"], "%s_transcripts.stats.gtf" % (chrom))        
        taskinfo = (chrom, params, output_file)
        tasklist.append(taskinfo)
    # create a pool of workers to process the tasks
    pool = multiprocessing.Pool(processes=processors)
    # submit tasks to the pool
    result_iter = pool.imap(measure_transcripts_worker, tasklist)
    # prevent any more tasks from being submitted 
    pool.close()
    # process results
    for task_outfile in result_iter:
        # copy all results to single file
        shutil.copyfileobj(open(task_outfile), outfileh)        
        # delete individual chromosome files
        if not keep_tmp:
            os.remove(task_outfile)
    return 0

if __name__ == '__main__':
    from optparse import OptionParser
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s',
                        level=logging.DEBUG)
    optionparser = OptionParser("usage: %prog [options] <config.xml>")
    optionparser.add_option('-p', '--processors', type="int", default=1)
    optionparser.add_option('--keep-tmp', action="store_true", default=False)
    (options, args) = optionparser.parse_args()
    config_file = args[0]
    outfileh = sys.stdout
    # read relevant parts in config file
    measure_transcripts_parallel(config_file, outfileh, 
                                 options.processors,
                                 keep_tmp=True)
    #params = parse_config_file(config_file)    
    #measure_transcripts('chrY', params, outfileh)
