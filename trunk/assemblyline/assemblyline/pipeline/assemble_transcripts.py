'''
Created on Dec 2, 2011

@author: mkiyer
'''
import argparse
import logging
import os

from assemblyline.lib.transcript_parser import parse_gtf
from assemblyline.lib.transcript import Exon, strand_int_to_str, NEG_STRAND, POS_STRAND
from assemblyline.lib.transcript_graph import TranscriptGraph
from assemblyline.lib.gtf import GTFFeature

def write_bed(chrom, name, strand, score, exons):
    #print "EXONS TO PRINT", exons    
    assert all(exons[0].start < x.start for x in exons[1:])
    assert all(exons[-1].end > x.end for x in exons[:-1])
    tx_start = exons[0].start
    tx_end = exons[-1].end    
    block_sizes = []
    block_starts = []
    for e in exons:
        block_starts.append(e.start - tx_start)
        block_sizes.append(e.end - e.start)        
    # make bed fields
    fields = [chrom, 
              str(tx_start), 
              str(tx_end),
              str(name),
              str(score),
              strand_int_to_str(strand),
              str(tx_start),
              str(tx_start),
              '0',
              str(len(exons)),
              ','.join(map(str,block_sizes)) + ',',
              ','.join(map(str,block_starts)) + ',']
    return fields

def write_gtf(chrom, gene_id, tx_id, strand, score, exons):
    # write transcript feature
    f = GTFFeature()
    f.seqid = chrom
    f.source = "AssemblyLine"
    f.feature_type = "transcript"
    f.start = min(e[0] for e in exons)
    f.end = max(e[1] for e in exons)
    f.score = score
    f.strand = strand
    f.phase = "."
    attrs = {}
    attrs['gene_id'] = gene_id
    attrs['tx_id'] = tx_id
    f.attrs = attrs
    return f

def get_transcript_id_label_map(transcripts):
    '''
    returns a dictionary with transcript ID keys and
    tuples of (transcript, exon_number) values
    '''
    id_map = {}
    for t in transcripts:
        # add scores to the score lookup table
        id_map[t.id] = t.label
    return id_map

def collapse_contiguous_nodes(path, strand):
    newpath = []    
    n = Exon(path[0].start, path[0].end)
    for i in xrange(1, len(path)):
        if strand == NEG_STRAND and (n.start == path[i].end):
            n.start = path[i].start
        elif n.end == path[i].start:
            n.end = path[i].end
        else:
            newpath.append(n)
            n = Exon(path[i].start, path[i].end)
    newpath.append(n)
    return newpath

def run(gtf_file, overhang_threshold, fraction_major_isoform, max_paths):
    locus_id = 1
    gene_id = 1
    tss_id = 1
    tx_id = 1
    for locus_transcripts in parse_gtf(open(gtf_file)):
        # gather properties of locus
        locus_start = locus_transcripts[0].start
        locus_end = max(tx.end for tx in locus_transcripts)
        logging.debug("Locus with %d transcripts chrom=%s start=%d end=%d" % 
                      (len(locus_transcripts), locus_transcripts[0].chrom,
                       locus_transcripts[0].start, locus_transcripts[-1].end))
        # build and refine transcript graph        
        transcript_graph = \
            TranscriptGraph.from_transcripts(locus_transcripts, 
                                             overhang_threshold)
        continue
        #transcript_graph.collapse(overhang_threshold=overhang_threshold)
        chrom = locus_transcripts[0].chrom
        # get transcript id -> transcript object mappings
        tx_id_label_map = get_transcript_id_label_map(locus_transcripts)
        # find isoforms of transcript
        for strand, loc_gene_id, loc_tss_id, score, path in \
            transcript_graph.assemble(max_paths=max_paths,
                                      fraction_major_path=fraction_major_isoform):
            # get transcripts and their sample labels from each path
            tx_ids = set().union(*[transcript_graph.get_exon_ids(n) for n in path])
            tx_labels = set(tx_id_label_map[id] for id in tx_ids)
            # use locus/gene/tss/transcript id to make gene name
            gene_name = "L%07d|G%07d|TSS%07d|TU%07d" % (locus_id,
                                                        gene_id + loc_gene_id, 
                                                        tss_id + loc_tss_id,
                                                        tx_id)
            # collapse contiguous nodes
            path = collapse_contiguous_nodes(path, strand) 
            # fix path
            if strand == NEG_STRAND:
                path.reverse()           
            fields = write_bed(chrom, gene_name, strand, score, path)
            # append extra columns to BED format
            #fields.append(str(len(tx_ids)))
            #fields.append(str(len(tx_labels)))
            print '\t'.join(fields)
            tx_id += 1
        gene_id += loc_gene_id + 1
        tss_id += loc_tss_id + 1 
        locus_id += 1

def main():
    # parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("-v", action="store_true", dest="verbose", default=False)
    parser.add_argument("--overhang", type=int, dest="overhang_threshold", default=100)
    parser.add_argument("--fraction-major-isoform", dest="fraction_major_isoform", type=float, default=0.10)
    parser.add_argument("--max-paths", dest="max_paths", type=int, default=1000)
    parser.add_argument("--score-attr", dest="score_attr", default="FPKM")
    parser.add_argument("filename")
    args = parser.parse_args()

    # set logging level
    if args.verbose:
        level = logging.DEBUG
    else:
        level = logging.WARNING
    logging.basicConfig(level=level,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")    

    # start algorithm
    run(args.filename, 
        args.overhang_threshold,
        args.fraction_major_isoform,
        args.max_paths)

if __name__ == '__main__': main()