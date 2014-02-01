'''
Created on Jan 8, 2014

@author: mkiyer
'''
import os
import sys
import argparse
import logging
import subprocess
import itertools
import shutil
from multiprocessing import Process, Queue

import pysam

from assemblyline.lib.base import which
from assemblyline.lib.seq import DNA_reverse_complement
from assemblyline.lib.transcript import parse_gtf, strand_int_to_str, strand_str_to_int, NEG_STRAND, NO_STRAND, POS_STRAND, Exon
from assemblyline.lib.batch_sort import batch_sort

SIGNALP_HEADER = ['name', 'Cmax', 'Cpos', 'Ymax', 'Ypos', 'Smax', 'Spos', 
                  'Smean', 'D', '?', 'Dmaxcut', 'Networksused']
PFAM_HEADER = ['seq_id', 'alignment_start', 'alignment_end', 
               'envelope_start', 'envelope_end', 'hmm_acc', 'hmm_name', 
               'type', 'hmm_start', 'hmm_end', 'hmm_length', 'bit_score',
               'e_value', 'significance', 'clan']

# batch sort configuration
SORT_BUFFER_SIZE = 32000

PFAM_FILES = ['Pfam-A.hmm',
              'Pfam-A.hmm.dat',
              'Pfam-A.hmm.h3f',
              'Pfam-A.hmm.h3i',
              'Pfam-A.hmm.h3m',
              'Pfam-A.hmm.h3p',
              'Pfam-B.hmm',
              'Pfam-B.hmm.dat',
              'Pfam-B.hmm.h3f',
              'Pfam-B.hmm.h3i',
              'Pfam-B.hmm.h3m',
              'Pfam-B.hmm.h3p']
def check_pfam_dir(path):
    for pfam_file in PFAM_FILES:
        if not os.path.exists(os.path.join(path, pfam_file)):
            return False
    return True

def parse_fasta_file(fasta_file):
    with open(fasta_file) as f:
        while True:
            lines = list(itertools.islice(f,2))
            if lines:
                yield lines
            else:
                break

def split_fasta_file(fasta_file, n, prefix):
    # open files
    filehandles = []
    for x in xrange(n):
        filename = '%s%d.fasta' % (prefix, x)
        filehandles.append(open(filename))
    x = 0
    for lines in parse_fasta_file(fasta_file):
        filehandles[x].writelines(lines)
        x  = (x + 1) % n
    # close files
    for f in filehandles:
        f.close()

GENETIC_CODE = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
    'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
}

class ORFInfo(object):
    FIELDS = ['transcript_id', 'gene_id', 'orf_id', 'frame', 
              'chrom', 'start', 'end', 'strand', 'exons', 'seq']
    SEQ_COL_NUM = 11
    
    @staticmethod
    def make_orf_id(x):
        return 'ORF%09d' % (x)
    @staticmethod
    def cmp_orf_id(x, y):
        x = int(x.split('ORF')[1])
        y = int(y.split('ORF')[1])
        return cmp(x, y)
    
    def __init__(self):
        for field in ORFInfo.FIELDS:
            setattr(self, field, None)
    
    def to_table(self):
        return [self.transcript_id,
                self.gene_id,
                str(self.orf_id),
                str(self.frame),
                str(len(self.seq)),
                self.chrom,
                str(self.start),
                str(self.end),
                self.strand,
                ','.join(map(str, (e.start for e in self.exons))),
                ','.join(map(str, (e.end for e in self.exons))),
                self.seq]

    @staticmethod
    def from_table(line):
        fields = line.strip().split('\t')
        self = ORFInfo()
        self.transcript_id = fields[0]
        self.gene_id = fields[1]
        self.orf_id = fields[2]
        self.frame = int(fields[3])
        self.chrom = fields[5]
        self.start = int(fields[6])
        self.end = int(fields[7])
        self.strand = fields[8]
        exon_starts = map(int, fields[9].split(','))
        exon_ends = map(int, fields[10].split(','))
        self.exons = [Exon(x,y) for x,y in zip(exon_starts, exon_ends)]
        self.seq = fields[ORFInfo.SEQ_COL_NUM]
        return self
    
    def to_bed(self, name=None):
        assert all(self.exons[0].start < x.start for x in self.exons[1:])
        assert all(self.exons[-1].end > x.end for x in self.exons[:-1])
        assert self.start == self.exons[0].start
        assert self.end == self.exons[-1].end
        block_sizes = []
        block_starts = []
        for e in self.exons:
            block_starts.append(e.start - self.start)
            block_sizes.append(e.end - e.start)
        if name is None:
            name = '%s|%s|frame=%d|length=%d' % (self.transcript_id, self.gene_id, self.frame, len(self.seq))
        # make bed fields
        fields = [self.chrom, 
                  str(self.start), 
                  str(self.end),
                  name,
                  '0',
                  self.strand,
                  str(self.start),
                  str(self.end),
                  '0',
                  str(len(self.exons)),
                  ','.join(map(str,block_sizes)) + ',',
                  ','.join(map(str,block_starts)) + ',']
        return fields

def get_transcript_dna_sequence(t, ref_fa):
    exon_seqs = []
    error_occurred = False
    for e in t.exons:
        seq = ref_fa.fetch(t.chrom, e.start, e.end)
        if (not seq) or (len(seq) < (e.end - e.start)):
            logging.warning("transcript %s exon %s:%d-%d not found in reference" % 
                            (t.attrs['transcript_id'], t.chrom, e.start, e.end))
            error_occurred = True
            break
        exon_seqs.append(seq)
    if error_occurred:
        return None
    # make fasta record
    seq = ''.join(exon_seqs).upper()
    # reverse complement negative stranded sequences
    if t.strand == NEG_STRAND:
        seq = DNA_reverse_complement(seq)
    return seq

def translate_dna_3frames(seq):
    proteins = [[] for x in xrange(3)]
    for n in range(0,len(seq),3):
        for frame in xrange(3):
            codon = seq[n+frame:n+frame+3]
            proteins[frame].append(GENETIC_CODE.get(codon, '.'))
    for frame in xrange(3):
        proteins[frame] = ''.join(proteins[frame]).rstrip('.')
    return proteins

def find_orfs(seq, start_codon='M', stop_codon='*'):
    '''
    start is inclusive and end is exclusive
    includes stop codon
    '''
    # find first start codon
    start = seq.find(start_codon)
    while (start >= 0):
        end = seq.find(stop_codon, start)
        if end == -1:
            yield start, len(seq), seq[start:]
            return
        yield start, end+1, seq[start:end+1]
        # find next start codon
        start = seq.find(start_codon, end+1)

def transcript_to_genome_pos(t, pos, reverse=False):
    reverse = (t.strand == NEG_STRAND) or reverse
    if reverse:
        exon_iter = xrange(len(t.exons)-1, -1, -1)
    else:
        exon_iter = xrange(0, len(t.exons))
    offset = 0
    for i in exon_iter:
        exon = t.exons[i]
        exon_size = exon.end - exon.start
        if pos < offset + exon_size:            
            if reverse:
                return exon.start + exon_size - (pos - offset) - 1
            else:
                return exon.start + (pos - offset)
        offset += exon_size
    if reverse:
        return t.start
    else:
        return (t.end-1)

def genome_interval_to_exons(start, end, t):
    assert start >= t.start
    assert end <= t.end    
    assert start < end
    i = 0
    while t.exons[i].end < start:
        i += 1
    start_exon = i
    while t.exons[i].end < end:
        i += 1
    end_exon = i
    newexons = []
    if start_exon == end_exon:
        newexons.append(Exon(start,end))
    else:
        newexons.append(Exon(start, t.exons[start_exon].end))
        for i in xrange(start_exon+1, end_exon):
            newexons.append(t.exons[i])
        newexons.append(Exon(t.exons[end_exon].start, end))
    return newexons

def orf_to_genome(t, frame, aa_start, aa_end):
    t_orf_start = frame + (3 * aa_start)
    t_orf_end = frame + (3 * aa_end) - 1 # make end position inclusive
    reverse = ((t.strand == NEG_STRAND) or (frame >= 3))
    g_orf_start = transcript_to_genome_pos(t, t_orf_start, reverse)
    g_orf_end = transcript_to_genome_pos(t, t_orf_end, reverse)
    if reverse:
        g_orf_start, g_orf_end = g_orf_end, g_orf_start    
    g_orf_end += 1 # end position is inclusive so adjust
    orf_exons = genome_interval_to_exons(g_orf_start, g_orf_end, t)
    return g_orf_start, g_orf_end, orf_exons

def translate_transcript(t, ref_fa):
    # get DNA sequence of transcript
    seq = get_transcript_dna_sequence(t, ref_fa)
    # translate on all reading frames
    if t.strand == NO_STRAND:
        aa_seqs = translate_dna_3frames(seq)
        aa_seqs.extend(translate_dna_3frames(DNA_reverse_complement(seq)))
    else:
        aa_seqs = translate_dna_3frames(seq)
    return aa_seqs

def to_fasta(name, seq):
    return '>%s\n%s' % (name, seq)

def get_transcript_orfs(t, aa_seqs, min_orf_length):
    for frame, aa_seq in enumerate(aa_seqs):
        for aa_start, aa_end, orf in find_orfs(aa_seq):
            if len(orf) < min_orf_length:
                continue
            orf_start, orf_end, orf_exons = \
                orf_to_genome(t, frame, aa_start, aa_end)
            orfinfo = ORFInfo()
            orfinfo.transcript_id = t.attrs['transcript_id']
            orfinfo.gene_id = t.attrs['gene_id']
            orfinfo.frame = frame
            orfinfo.chrom = t.chrom
            orfinfo.start = orf_start
            orfinfo.end = orf_end
            orfinfo.strand = strand_int_to_str(t.strand)
            orfinfo.exons = orf_exons
            orfinfo.seq = orf
            yield orfinfo

def group_unique_orfs(line_iter):
    cur_orfs = []
    num_unique_orfs = 1
    cur_orf_id = ORFInfo.make_orf_id(num_unique_orfs)
    for line in line_iter:
        orf = ORFInfo.from_table(line)
        if (len(cur_orfs) > 0) and (cur_orfs[0].seq != orf.seq):
            yield cur_orfs
            cur_orfs = []
            num_unique_orfs += 1
            cur_orf_id = ORFInfo.make_orf_id(num_unique_orfs)
        orf.orf_id = cur_orf_id
        cur_orfs.append(orf)
    if len(cur_orfs) > 0:
        yield cur_orfs

def run_signalp(fasta_files, output_file, tmp_dir):
    def _consumer(worker_index, worker_fasta_file, output_queue):
        # process with signalp
        logging.debug('[Worker %d] Running SignalP' % (worker_index))
        args = ['signalp', '-f', 'short', '-t', 'euk', worker_fasta_file]
        tmp_output_file = os.path.join(tmp_dir, 'worker%d.signalp.txt' % (worker_index))
        with open(tmp_output_file, 'w') as outfileh:
            retcode = subprocess.call(args, stdout=outfileh)
            if retcode != 0:
                logging.error('[Worker %d] Error running SignalP' % (worker_index))
        # put results into output queue
        with open(tmp_output_file) as f:
            for line in f:
                if not line:
                    continue
                line = line.strip()
                if not line:
                    continue
                if line.startswith('#'):
                    continue
                output_queue.put(line)
        output_queue.put(None)
        logging.debug('[Worker %d] Finished' % (worker_index))
    # create multiprocessing queue for passing data
    output_queue = Queue(maxsize=len(fasta_files)*4) 
    # start consumer processes
    procs = []
    for i in xrange(len(fasta_files)):
        p = Process(target=_consumer, args=(i, fasta_files[i], output_queue))
        p.start()
        procs.append(p)
    # get results from consumers
    num_alive = len(procs)
    tmp_output_file = os.path.join(tmp_dir, 'signalp.merge.txt')
    with open(tmp_output_file, 'w') as f:
        while num_alive > 0:
            result = output_queue.get()
            if result is None:
                num_alive -= 1
                logging.debug("Main process detected worker finished, %d still alive" % (num_alive))
            else:
                print >>f, result
    logging.debug("Joining all processes")
    # wait for consumers to finish
    for p in procs:
        p.join()
    # sort
    logging.debug('Sorting signal peptide results')
    def sort_signalp(line):
        return line.split(None,1)[0]
    sort_tmp_dir = os.path.join(tmp_dir, 'sort_tmp')
    os.makedirs(sort_tmp_dir)
    batch_sort(input=tmp_output_file,
               output=output_file,
               key=sort_signalp,
               buffer_size=SORT_BUFFER_SIZE,
               tempdirs=[sort_tmp_dir])
    # cleanup temporary files
    shutil.rmtree(sort_tmp_dir)
    return 0

def run_pfam(fasta_files, pfam_dir, output_file, tmp_dir):
    def _consumer(worker_index, worker_fasta_file, output_queue):
        # process with pfam
        logging.debug('[Worker %d] Running pfam_scan.pl' % (worker_index))
        tmp_output_file = os.path.join(tmp_dir, 'worker%d.pfam.txt' % (worker_index))
        args = ['pfam_scan.pl', 
                '-cpu', '1',
                '-pfamB', 
                '-fasta', worker_fasta_file,
                '-dir', pfam_dir,
                '-outfile', tmp_output_file]
        retcode = subprocess.call(args)
        if retcode != 0:
            logging.error('[Worker %d] Error running pfam_scan.pl' % (worker_index))
        # put results into output queue
        with open(tmp_output_file) as f:
            for line in f:
                if not line:
                    continue
                line = line.strip()
                if not line:
                    continue
                if line.startswith('#'):
                    continue
                output_queue.put(line)
        output_queue.put(None)
        logging.debug('[Worker %d] Finished' % (worker_index))
    # create multiprocessing queues for passing data
    output_queue = Queue(maxsize=len(fasta_files)*4) 
    # start consumer processes
    procs = []
    for i in xrange(len(fasta_files)):
        p = Process(target=_consumer, args=(i, fasta_files[i], output_queue))
        p.start()
        procs.append(p)
    # get results from consumers
    num_alive = len(procs)
    tmp_output_file = os.path.join(tmp_dir, 'pfam.merge.txt')
    with open(tmp_output_file, 'w') as f:
        while num_alive > 0:
            result = output_queue.get()
            if result is None:
                num_alive -= 1
                logging.debug("Main process detected worker finished, %d still alive" % (num_alive))
            else:
                print >>f, result
    logging.debug("Joining all processes")
    # wait for consumers to finish
    for p in procs:
        p.join()
    # sort
    logging.debug('Sorting pfam results')
    def sort_pfam(line):
        return line.split(None,1)[0]
    sort_tmp_dir = os.path.join(tmp_dir, 'sort_tmp')
    os.makedirs(sort_tmp_dir)
    batch_sort(input=tmp_output_file,
               output=output_file,
               key=sort_pfam,
               buffer_size=SORT_BUFFER_SIZE,
               tempdirs=[sort_tmp_dir])
    # cleanup temporary files
    shutil.rmtree(sort_tmp_dir)
    return 0

def merge_results(orf_table_file, signalp_file, pfam_file, output_file):
    class SortedFileParser:
        def __init__(self, fileh, key_col=0, key_cmp=ORFInfo.cmp_orf_id, sep=None):
            self.fileh = fileh
            self.key_col = key_col
            self.key_cmp = key_cmp
            self.fields = None
            self.sep = sep
        def seek(self, orf_id):
            try:
                while True:
                    if self.fields is None:
                        c = -1
                    else:
                        c = self.key_cmp(self.fields[self.key_col], orf_id)
                    if c == 0:
                        # found matching line
                        return self.fields
                    elif c == 1:
                        # larger than current id so no match
                        return None
                    # smaller than current id so keep seeking
                    self.fields = self.fileh.next().strip().split(self.sep)
            except StopIteration:
                pass
    with open(output_file, 'w') as outfile:
        with open(orf_table_file) as infile:
            pfam_parser = SortedFileParser(open(pfam_file), key_col=0)
            signalp_parser = SortedFileParser(open(signalp_file), key_col=0)
            for line in infile:
                orfinfo = ORFInfo.from_table(line)
                signalp_fields = signalp_parser.seek(orfinfo.orf_id)
                if signalp_fields is None:
                    signalp_fields = ['None'] * len(SIGNALP_HEADER)
                pfam_fields = pfam_parser.seek(orfinfo.orf_id)
                if pfam_fields is None:
                    pfam_fields = ['None'] * len(PFAM_HEADER)
                fields = orfinfo.to_table()
                fields.extend(signalp_fields)
                fields.extend(pfam_fields)
                print >>outfile, '\t'.join(fields)

def orf_analysis(gtf_file, genome_fasta_file, pfam_dir, 
                 output_dir, min_orf_length, num_processes):
    #
    # extract transcript DNA sequences, translate to protein, and
    # search for ORFs
    #
    logging.debug('Finding ORFs in transcript sequences')
    # output files
    tmp_dir = os.path.join(output_dir, 'tmp')
    if not os.path.exists(tmp_dir):
        logging.info("Creating tmp directory '%s'" % (tmp_dir))
        os.makedirs(tmp_dir)
    orf_bed_file = os.path.join(output_dir, 'transcript_orfs.bed')
    unique_orf_file = os.path.join(output_dir, 'unique_orfs.txt')
    unique_orf_bed_file = os.path.join(output_dir, 'unique_orfs.bed')
    orf_file = os.path.join(tmp_dir, 'transcript_orfs.no_ids.txt')
    sorted_orf_file = os.path.join(tmp_dir, 'transcript_orfs.no_ids.sortbyorf.txt')
    sorted_orf_id_file = os.path.join(tmp_dir, 'transcript_orfs.sortbyorf.txt')
    signalp_file = os.path.join(output_dir, 'signalp.txt')
    pfam_file = os.path.join(output_dir, 'pfam.txt')
    merged_orf_id_file = os.path.join(output_dir, 'transcript_orfs.sortbyorf.merged.txt')
    sorted_merged_orf_id_file = os.path.join(output_dir, 'transcript_orfs.sortbytranscript.merged.txt')
    # open output files
    orf_fileh = open(orf_file, 'w')
    orf_bed_fileh = open(orf_bed_file, 'w')
    # open genome fasta file
    ref_fa = pysam.Fastafile(genome_fasta_file)
    num_finished = 1
    for locus_transcripts in parse_gtf(open(gtf_file)):
        for t in locus_transcripts:
            # get amino acid sequences in all reading frames
            aa_seqs = translate_transcript(t, ref_fa)
            # get ORFs
            for orfinfo in get_transcript_orfs(t, aa_seqs, min_orf_length):
                print >>orf_fileh, '\t'.join(orfinfo.to_table())
                print >>orf_bed_fileh, '\t'.join(orfinfo.to_bed())
            if (num_finished % 10000) == 0:
                logging.debug('Processed %d transcripts' % (num_finished))
            num_finished += 1
    # cleanup
    orf_fileh.close()
    orf_bed_fileh.close()
    #
    # sort ORF table by ORF amino acid sequence to group identical ORFs 
    # together
    #
    logging.debug('Sorting ORFs by amino acid sequence')
    def sort_by_seq(line):
        '''comparison function for batch_sort'''
        fields = line.strip().split('\t')
        return fields[ORFInfo.SEQ_COL_NUM]
    sort_tmp_dir = os.path.join(tmp_dir, 'sort_tmp')
    os.makedirs(sort_tmp_dir)
    batch_sort(input=orf_file,
               output=sorted_orf_file,
               key=sort_by_seq,
               buffer_size=SORT_BUFFER_SIZE,
               tempdirs=[sort_tmp_dir])
    shutil.rmtree(sort_tmp_dir)
    #
    # assign each ORF a unique id and write to FASTA file
    #
    logging.debug('Determining unique ORFs')
    orf_fasta_prefix = os.path.join(tmp_dir, 'orf')
    orf_fasta_files = []
    orf_fasta_sizes = []
    for i in xrange(num_processes):
        orf_fasta_files.append(open('%s%d.fasta' % (orf_fasta_prefix, i), 'w'))
        orf_fasta_sizes.append(0)
    orf_file_index = 0
    outfileh = open(sorted_orf_id_file, 'w')
    unique_orf_fileh = open(unique_orf_file, 'w')
    print >>unique_orf_fileh, '\t'.join(['orf_id', 'orf_length', 'total_occurrences', 'unique_genomic_occurrences'])
    unique_orf_bed_fileh = open(unique_orf_bed_file, 'w')
    with open(sorted_orf_file) as infileh:
        for orfs in group_unique_orfs(infileh):
            # write to master transcript/ORF table
            for orf in orfs:
                print >>outfileh, '\t'.join(orfinfo.to_table())
            # write ORF to fasta file
            lines = to_fasta(orfs[0].orf_id, orfs[0].seq.strip('*'))
            print >>orf_fasta_files[orf_file_index], lines
            orf_fasta_sizes[orf_file_index] += 1
            # advance to next fasta file
            orf_file_index = (orf_file_index + 1) % (num_processes)
            # group by genomic position and write ORFs to BED file
            unique_genome_orfs = {}
            for orf in orfs:
                k = (orf.chrom, orf.strand, tuple(orf.exons))
                if k in unique_genome_orfs:
                    continue
                unique_genome_orfs[k] = orf
            for orf in unique_genome_orfs.itervalues():
                print >>unique_orf_bed_fileh, '\t'.join(orf.to_bed(orf.orf_id))
            # write unique ORF to tab-delimited text file
            fields = [orfs[0].orf_id, len(orfs[0].seq), len(orfs), len(unique_genome_orfs)]
            print >>unique_orf_fileh, '\t'.join(map(str,fields))
    # cleanup
    unique_orf_bed_fileh.close()
    outfileh.close()
    # get fasta files with lines written
    orf_fasta_file_names = []
    for i in xrange(len(orf_fasta_files)):
        orf_fasta_files[i].close()
        if orf_fasta_sizes[i] > 0:
            orf_fasta_file_names.append(orf_fasta_files[i].name)
    #
    # search FASTA file against signalp
    #
    logging.debug('Searching for signal peptides')
    retcode = run_signalp(orf_fasta_file_names, signalp_file, tmp_dir)
    if retcode != 0:
        logging.error('Error searching for signal peptides')
        return 1 
    #
    # search FASTA file against Pfam
    #
    logging.error('Scanning for Pfam domains')
    retcode = run_pfam(orf_fasta_file_names, pfam_dir, pfam_file, tmp_dir)
    if retcode != 0:
        logging.error('Error running pfam_scan.pl')
    #
    # merge results from Pfam and signalp
    #
    logging.debug('Merging SignalP and Pfam results')
    merge_results(sorted_orf_id_file, signalp_file, pfam_file, merged_orf_id_file)
    #
    # sort by transcript id
    #
    logging.debug('Sorting ORFs by transcript ID')
    def sort_by_transcript_id(line):
        return line.split('\t', 1)[0]
    sort_tmp_dir = os.path.join(tmp_dir, 'sort_tmp')
    os.makedirs(sort_tmp_dir)
    batch_sort(input=merged_orf_id_file,
               output=sorted_merged_orf_id_file,
               key=sort_by_transcript_id,
               buffer_size=SORT_BUFFER_SIZE,
               tempdirs=[sort_tmp_dir])
    shutil.rmtree(sort_tmp_dir)
    # cleanup
    ref_fa.close()
    return 0

def full_transcript_analysis(gtf_file, genome_fasta_file, pfam_dir, 
                             output_dir, num_processes):
    # output files    
    pfam_file = os.path.join(output_dir, 'full_length_pfam.txt')
    # open genome fasta file
    ref_fa = pysam.Fastafile(genome_fasta_file)
    # convert transcripts to amino acid sequences and write to fasta file
    logging.debug('Writing transcript amino acid sequences to FASTA file(s)')
    tmp_dir = os.path.join(output_dir, 'tmp')
    os.makedirs(tmp_dir)
    fasta_prefix = os.path.join(tmp_dir, 'full')
    fasta_files = []
    fasta_sizes = []
    for i in xrange(num_processes):
        fasta_files.append(open('%s%d.fasta' % (fasta_prefix, i), 'w'))
        fasta_sizes.append(0)
    num_finished = 1
    fasta_file_index = 0
    for locus_transcripts in parse_gtf(open(gtf_file)):
        for t in locus_transcripts:
            # get amino acid sequences in all reading frames
            aa_seqs = translate_transcript(t, ref_fa)
            for frame, aa_seq in enumerate(aa_seqs):
                lines = to_fasta('%s|frame=%d' % (t.attrs['transcript_id'], frame), aa_seq)
                print >>fasta_files[fasta_file_index], lines
            fasta_sizes[fasta_file_index] += 1
            fasta_file_index = (fasta_file_index + 1) % num_processes
            if (num_finished % 10000) == 0:
                logging.debug('Processed %d transcripts' % (num_finished))
            num_finished += 1
    # get fasta files with lines written
    fasta_file_names = []
    for i in xrange(len(fasta_files)):
        fasta_files[i].close()
        if fasta_sizes[i] > 0:
            fasta_file_names.append(fasta_files[i].name)  
    # cleanup
    ref_fa.close()
    #
    # search FASTA file against Pfam
    #
    logging.debug('Scanning for Pfam domains')
    retcode = run_pfam(fasta_file_names, pfam_dir, pfam_file, tmp_dir)
    if retcode != 0:
        logging.error('Error running pfam_scan.pl')
        return retcode
    return 0

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument('--pfam', dest='pfam_dir')
    parser.add_argument('--genome-fasta', dest='genome_fasta_file')
    parser.add_argument('--min-orf-length', dest='min_orf_length', type=int, default=30)
    parser.add_argument('-o', '--output-dir', dest='output_dir', default='out')
    parser.add_argument('-p', '--num-processes', dest='num_processes', type=int, default=1)
    parser.add_argument('--mode', dest='mode', choices=['orf', 'full'], default='orf')
    parser.add_argument('gtf_file')
    args = parser.parse_args()
    # get args
    pfam_dir = args.pfam_dir
    genome_fasta_file = args.genome_fasta_file
    gtf_file = args.gtf_file
    min_orf_length = args.min_orf_length
    num_processes = args.num_processes
    output_dir = args.output_dir
    mode = args.mode
    # check command line parameters
    if which('pfam_scan.pl') is None:
        parser.error("'pfam_scan.pl' not found in PATH")
    if which('hmmscan') is None:
        parser.error("'hmmscan' not found in PATH")
    if which('signalp') is None:
        parser.error("'signalp' not found in PATH")
    if not check_pfam_dir(pfam_dir):
        parser.error("Required Pfam-A and Pfam-B files not found at '%s'" % (pfam_dir))
    if not os.path.exists(genome_fasta_file):
        parser.error("Genome FASTA file '%s' not found" % (genome_fasta_file))
    if not os.path.exists(gtf_file):
        parser.error("GTF file '%s' not found" % (gtf_file))
    #if os.path.exists(output_dir):
    #    parser.error("Output directory '%s' already exists" % (output_dir))
    # create output dir
    if not os.path.exists(output_dir):
        logging.info("Creating output directory '%s'" % (output_dir))
        os.makedirs(output_dir)
    if mode == 'orf':    
        return orf_analysis(gtf_file, genome_fasta_file, pfam_dir, 
                            output_dir, min_orf_length, num_processes)
    else:
        return full_transcript_analysis(gtf_file, genome_fasta_file, 
                                        pfam_dir, output_dir, num_processes)
    return 0


if __name__ == '__main__':
    sys.exit(main())