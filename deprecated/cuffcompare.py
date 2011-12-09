'''
Created on Feb 17, 2010

@author: mkiyer
'''

from optparse import OptionParser
import os
import glob
import logging
import collections
import numpy as np
import multiprocessing

import genome
import gtf2

import pipelineconfig
from transcript import Transcript

cuffcompare_combined_attrs = {'gene_id': str,
                              'transcript_id': str,
                              'exon_number': int,
                              'oId': str,
                              'nearest_ref': str,
                              'class_code': str,
                              'tss_id': str}

stats_gtf_attrs = dict(cuffcompare_combined_attrs.items())
stats_gtf_own_attrs = [('ann', int),
                       ('max_ann', int),
                       ('recur', int),
                       ('fpkm_lo', float),
                       ('fpkm_hi', float),
                       ('fpkm', float),
                       ('uavg', float),
                       ('ustd', float),
                       ('ulong', int),
                       ('max_recur', int),
                       ('max_fpkm_lo', float),
                       ('max_fpkm_hi', float),
                       ('max_fpkm', float),
                       ('max_length', int)]
stats_gtf_attrs.update(stats_gtf_own_attrs)

classify_gtf_attrs = dict(stats_gtf_attrs.items())
classify_gtf_attrs.update({'pred': int})


def parse_gtf_transcript_features(fileh, attrs=None):
    # read first feature to initialize parser
    feature_iter = gtf2.read_gtf2_file(fileh, attrs)
    feature = feature_iter.next()
    transcript_features = [feature]
    # keep track of current gene/transcript
    current_transcript_id = feature.attrs['transcript_id']    
    # now that initialized using first line read rest of 
    # transcripts
    for feature in feature_iter:
        # get gene/transcript ids
        #next_gene_id = feature.attrs['gene_id']
        transcript_id = feature.attrs['transcript_id']
        if current_transcript_id != transcript_id:
            # build transcript from current set of features
            yield transcript_features
            # initialize new transcript
            current_transcript_id = transcript_id
            transcript_features = []
        # add features to current transcript
        transcript_features.append(feature)
    # clean up last transcript in data
    if len(transcript_features) > 0:
        yield transcript_features

def parse_combined_gtf_features(fileh):
    return parse_gtf_transcript_features(fileh, cuffcompare_combined_attrs)

def parse_combined_gtf_transcripts(fileh):
    for transcript_features in parse_combined_gtf_features(fileh):
        transcript = Transcript(transcript_features[0].attrs['transcript_id'], 
                                transcript_features[0].seqid,
                                transcript_features[0].strand,
                                [(f.start, f.end) for f in transcript_features])
        yield transcript

class TrackingData(object):
    pass

def parse_tracking_file(fileh, ordered_sample_names):
    for line in fileh:
        fields = line.strip().split('\t')
        # extract transcript information
        tx_id = fields[0]
        locus_id = fields[1]
        ref_ids = fields[2].strip()
        if ref_ids != '-':
            x, y = ref_ids.split('|')
            ref_gene_id = x
            ref_tx_id = y
        else:
            ref_gene_id = None
            ref_tx_id = None
        class_code = fields[3]
        num_inputs = len(fields[4:])
        # make arrays to store sample information
        fmi_arr = 100 * np.ones(num_inputs)
        fpkm_arr = np.zeros(num_inputs)
        fpkm_lo_arr = np.zeros(num_inputs)
        fpkm_hi_arr = np.zeros(num_inputs)
        cov_arr = np.zeros(num_inputs)
        #recurrence = 0
        # extract information about individual input samples
        #for field in fields[4:]:
        samples_found = set()
        for index, field in enumerate(fields[4:]):
            # qJ:<gene_id>|<transcript_id>|<FMI>|<RPKM>|<conf_lo>|<conf_hi>|<cov>
            if field == '-':
                continue
            index_string, sample_info = field.split(':')
            #index = (int(index_string[1:]) - 1)
            assert index == (int(index_string[1:]) - 1)
            # measure recurrence on technical replicate level
            #recurrence += 1
            # measure recurrence on a biological replicate level
            sample_name = ordered_sample_names[index]
            samples_found.add(sample_name)
            # measure properties of expression and coverage            
            sample_fields = sample_info.split('|')
            sample_gene_id = sample_fields[0]
            sample_tx_id = sample_fields[1]            
            fmi_arr[index] = float(sample_fields[2])
            fpkm_arr[index] = float(sample_fields[3])
            fpkm_lo_arr[index] = float(sample_fields[4])
            fpkm_hi_arr[index] = float(sample_fields[5])
            cov_arr[index] = float(sample_fields[6])
        # return info
        tobj = TrackingData()
        tobj.transcript_id = tx_id
        tobj.locus_id = locus_id
        #tobj.recurrence = recurrence
        tobj.recurrence = len(samples_found)
        tobj.fmi = fmi_arr
        tobj.fpkm = fpkm_arr
        tobj.fpkm_lo = fpkm_lo_arr
        tobj.fpkm_hi = fpkm_hi_arr
        tobj.cov = cov_arr
        yield tobj

def extract_chroms_from_gtf(gtf_file, output_dir, suffix, genome_version):
    chroms = genome.get_genome(genome_version).get_chrom_names()
    chrom_file_names = {}
    chrom_fhds = {}
    for chrom in chroms:
        out_file = '_'.join([chrom, suffix]) + '.gtf'
        chrom_file_names[chrom] = out_file
        out_path = os.path.abspath(os.path.join(output_dir, out_file))
        chrom_fhds[chrom] = open(out_path, 'w')
    for line in open(gtf_file):
        chrom = line.strip().split('\t')[0]
        if chrom not in chroms:
            continue
        chrom_fhds[chrom].write(line)
    for chrom in chroms:
        chrom_fhds[chrom].close()
    return chrom_file_names

def cuffcompare_worker(task):
    chrom, cmdstring = task
    os.system(cmdstring)
    return chrom

def run_cuffcompare(config_file, processors=1):
    # read config
    config = pipelineconfig.PipelineConfig(config_file)
    logging.debug("Sample pools file: %s" % config.sample_pools_file)
    logging.debug("Tuxedo root dir: %s" % config.tuxedo_root_dir)        
    logging.debug("Cuffcompare path: %s" % config.cuffcompare_output_path)        
    logging.debug("Cuffcompare binary: %s" % config.cuffcompare_executable)        
    logging.debug("Cuffcompare reference GTF: %s" % config.cuffcompare_ref_gtf)

    # create output path
    if not os.path.exists(config.cuffcompare_output_path):
        os.makedirs(config.cuffcompare_output_path)

    # partition the GTF files by chromosome to compensate for memory problems 
    input_files = collections.defaultdict(lambda: [])
    for sample_name, library_ids in zip(config.sample_pools.names, config.sample_pools.library_ids):        
        for library_id in library_ids:
            path = os.path.join(config.tuxedo_data_dir, library_id)            
            if not os.path.exists(path):
                logging.error("%s %s analysis path %s does not exist" % (sample_name, library_id, path))
                continue
            gtf_file = os.path.join(path, 'transcripts.gtf')
            if not os.path.exists(gtf_file):
                logging.error("%s %s cufflinks gtf file does not exist" % (sample_name, library_id))
                continue
            logging.debug("%s: Partitioning by chromosome.." % (path))
            chrom_files = extract_chroms_from_gtf(gtf_file, config.cuffcompare_output_path, library_id, config.genome_version)
            for k,v in chrom_files.iteritems():
                input_files[k].append(v)

    # change directory to the output path
    curdir = os.getcwd()
    os.chdir(config.cuffcompare_output_path)
    # run cuffcompare once per chromosome
    tasklist = []
    for chrom in genome.get_genome(config.genome_version).get_chrom_names():
        logging.debug("Running cuffcompare on chromosome %s" % (chrom))
        # make file containing GTF files
        inp_file_name = chrom + '.input_files'
        f = open(inp_file_name, 'w')
        f.write('\n'.join(input_files[chrom]))
        f.close()
        cmd = [config.cuffcompare_executable, '-T', '-o', chrom + '.stats']
        if config.cuffcompare_ref_enabled:
            cmd.extend(['-r', config.cuffcompare_ref_gtf, '-R'])
        cmd.append('-i ' + inp_file_name)
        cmdstring = ' '.join(cmd)
        logging.debug('Cuffcompare command: ' + cmdstring)
        tasklist.append((chrom, cmdstring))

    # create a pool of workers to process the tasks
    pool = multiprocessing.Pool(processes=processors)
    # submit tasks to the pool
    result_iter = pool.imap(cuffcompare_worker, tasklist)
    # prevent any more tasks from being submitted 
    pool.close()
    # process results
    for chrom in result_iter:
        logging.debug("Finished cuffcompare on %s" % chrom)
        # remove chrom files
        for f in input_files[chrom]:
            os.remove(f)
    # change back to original directory
    os.chdir(curdir)
    return 0


if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)
    optionparser = OptionParser("usage: %prog [options] <config.xml>")
    optionparser.add_option('-p', '--processors', type="int", default=1)    
    (options, args) = optionparser.parse_args()
    config_file = args[0]
    run_cuffcompare(config_file, processors=options.processors)

