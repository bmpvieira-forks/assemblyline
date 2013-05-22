'''
Created on Oct 5, 2012

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
import collections
import sys
from multiprocessing import Process, JoinableQueue

# project imports
import assemblyline
import assemblyline.lib.config as config
from assemblyline.lib.bx.intersection import Interval, IntervalTree
from assemblyline.lib.gtf import parse_loci, merge_sort_gtf_files
from assemblyline.lib.transcript import transcripts_from_gtf_lines, \
    POS_STRAND, NEG_STRAND, NO_STRAND
from assemblyline.lib.base import Category, GTFAttr, FLOAT_PRECISION
from assemblyline.lib.assemble.transcript_graph import \
    find_exon_boundaries, split_exons

CInfo = collections.namedtuple('CategoryInfo',
                               ['category',
                                'ref',
                                'ann_intron_ratio',
                                'ann_cov_ratio',
                                'is_test'])

def compute_coverage_overlap(nodes1, nodes2):
    a = set(nodes1)
    b = set(nodes2)
    shared_nodes = a.intersection(b)
    union_nodes = a.union(b)
    shared_length = sum((n[1] - n[0]) for n in shared_nodes)    
    total_length = sum((n[1] - n[0]) for n in union_nodes)
    return shared_length, total_length

def find_best_coverage_overlap(nodes, reftuples, ignore_test=False):
    # find the reference transcript with the best overlap
    best_ref_t = None
    best_shared_cov_ratio = 0.0
    for ref_t, ref_nodes in reftuples:
        is_test = bool(int(ref_t.attrs[GTFAttr.TEST]))
        if is_test and ignore_test:
            continue
        shared_cov, union_cov = compute_coverage_overlap(nodes, ref_nodes)
        shared_ratio = float(shared_cov) / union_cov
        if shared_ratio > best_shared_cov_ratio:
            best_ref_t = ref_t
            best_shared_cov_ratio = shared_ratio
    return best_ref_t, best_shared_cov_ratio

def find_best_intron_overlap(test_nodes, test_introns, reftuples, 
                             ignore_test=False):
    best_ref_t = None
    best_shared_intron_ratio = 0.0
    best_shared_cov_ratio = 0.0    
    for ref_t, ref_nodes in reftuples:
        is_test = bool(int(ref_t.attrs[GTFAttr.TEST]))
        if is_test and ignore_test:
            continue
        ref_introns = set(ref_t.iterintrons())
        shared_introns = test_introns.intersection(ref_introns)
        union_introns = test_introns.union(ref_introns)
        shared_intron_ratio = float(len(shared_introns)) / len(union_introns)
        shared_cov, total_cov = compute_coverage_overlap(test_nodes, ref_nodes)
        shared_cov_ratio = float(shared_cov) / total_cov
        if ((shared_intron_ratio > best_shared_intron_ratio) or
            ((shared_intron_ratio == best_shared_intron_ratio) and
             (shared_cov_ratio > best_shared_cov_ratio))):
            best_ref_t = ref_t
            best_shared_intron_ratio = shared_intron_ratio
            best_shared_cov_ratio = shared_cov_ratio
    return best_ref_t, best_shared_intron_ratio, best_shared_cov_ratio

def categorize_transcript(t, nodes, introns, 
                          shared_intron_refs,
                          same_strand_refs,
                          opp_strand_refs,
                          intron_tree,
                          ignore_test=False):
    if len(shared_intron_refs) > 0:
        # find reference transcript with best intron overlap
        # and break ties using total coverage overlap
        best_ref_t, ann_intron_ratio, ann_cov_ratio = \
            find_best_intron_overlap(nodes, introns, 
                                     shared_intron_refs,                                         
                                     ignore_test)
        if best_ref_t is not None:    
            # determine whether this is a 'test' transcript
            is_test = bool(int(best_ref_t.attrs[GTFAttr.TEST]))
            return CInfo(category=Category.SAME_STRAND ,
                         ref=best_ref_t,
                         ann_cov_ratio=ann_cov_ratio,
                         ann_intron_ratio=ann_intron_ratio,
                         is_test=is_test)                
            
    if len(same_strand_refs) > 0:
        # find the reference transcript with the best overlap
        best_ref_t, ann_cov_ratio = \
            find_best_coverage_overlap(nodes, 
                                       same_strand_refs, 
                                       ignore_test)
        if best_ref_t is not None:
            # determine whether this is a 'test' transcript
            is_test = bool(int(best_ref_t.attrs[GTFAttr.TEST]))
            return CInfo(category=Category.SAME_STRAND,
                         ref=best_ref_t,
                         ann_cov_ratio=ann_cov_ratio,
                         ann_intron_ratio=0.0,
                         is_test=is_test)                

    # not a reference transcript
    best_ref_t = None
    is_test = False
    ann_cov_ratio = 0.0
    if len(opp_strand_refs) > 0:
        # transcript has coverage overlapping on the opposite strand
        # compared to reference transcripts
        category = Category.OPP_STRAND
        # find the reference transcript with the best overlap
        best_ref_t, ann_cov_ratio = \
            find_best_coverage_overlap(nodes, 
                                       opp_strand_refs, 
                                       ignore_test=False)
    else:
        # transcript has no coverage overlapping a reference transcript
        # so it must be either intronic, interleaving, or intergenic
        # search for introns overlapping transcript
        found_hit = False
        categories = set()
        for hit in intron_tree.find(t.start, t.end):
            if ((t.strand == hit.strand) and 
                ((hit.start,hit.end) in introns)):
                continue
            found_hit = True
            # check if there is an intron that encompasses the 
            # entire transcript
            if (hit.start < t.start) and (hit.end > t.end):
                if t.strand == NO_STRAND:
                    categories.add(Category.INTRONIC_AMBIGUOUS)
                    # no need to check other introns for unstranded 
                    break
                elif hit.strand == t.strand:
                    categories.add(Category.INTRONIC_SAME_STRAND)
                else:
                    categories.add(Category.INTRONIC_OPP_STRAND)
        if not found_hit:
            # no overlap with introns
            category = Category.INTERGENIC
        elif len(categories) == 1:
            category = categories.pop()
        elif len(categories) > 1:
            # overlaps introns on both strand
            category = Category.INTRONIC_AMBIGUOUS
        else:
            # a single intron does not encompass the transcript
            category = Category.INTERLEAVING
    return CInfo(category=category,
                 ref=best_ref_t,
                 ann_cov_ratio=ann_cov_ratio,
                 ann_intron_ratio=0.0,
                 is_test=is_test)

def compute_recurrence_and_score(nodes, node_data):
    # gather recurrence and scores across each node
    total_length = 0.0
    total_score = 0.0
    total_recur = 0.0
    total_pctrank = 0.0
    for n in nodes:
        nd = node_data[n]
        length = float(n[1] - n[0])
        total_score += nd['score'] * length
        total_pctrank += nd['pct'] * length
        total_recur += len(nd['ids']) * length       
        total_length += length
    # calculate statistics
    mean_score = total_score / total_length
    mean_pctrank = total_pctrank / total_length
    mean_recur = total_recur / total_length
    return mean_score, mean_pctrank, mean_recur

def resolve_strand(nodes, node_score_dict, ref_node_dict):
    # find strand with highest score
    total_scores = [0.0, 0.0]
    ref_bp = [0, 0]
    for n in nodes:
        length = (n[1] - n[0])
        if n in node_score_dict:
            scores = node_score_dict[n]
            total_scores[POS_STRAND] += (scores[POS_STRAND]*length)
            total_scores[NEG_STRAND] += (scores[NEG_STRAND]*length)
        if n in ref_node_dict:
            strand_ref_ids = ref_node_dict[n]
            if len(strand_ref_ids[POS_STRAND]) > 0:
                ref_bp[POS_STRAND] += length
            if len(strand_ref_ids[NEG_STRAND]) > 0:
                ref_bp[NEG_STRAND] += length
    if sum(total_scores) > FLOAT_PRECISION:
        if total_scores[POS_STRAND] >= total_scores[NEG_STRAND]:
            return POS_STRAND
        else:
            return NEG_STRAND
    if sum(ref_bp) > 0:
        if ref_bp[POS_STRAND] >= ref_bp[NEG_STRAND]:
            return POS_STRAND
        else:
            return NEG_STRAND
    return NO_STRAND

def annotate_locus(transcripts, 
                   gtf_sample_attr): 
    # store reference introns
    # (strand,start,end) -> ids (set) 
    ref_intron_dict = collections.defaultdict(lambda: [])
    ref_node_dict = collections.defaultdict(lambda: ([],[]))
    node_score_dict = collections.defaultdict(lambda: [0.0, 0.0])
    all_introns = set()
    # find the intron domains of the transcripts
    boundaries = find_exon_boundaries(transcripts)
    # add transcript to intron and graph data structures
    inp_transcripts = []
    for t in transcripts:
        # separate ref and nonref transcripts
        is_ref = bool(int(t.attrs[GTFAttr.REF]))
        if is_ref:
            # split exons that cross boundaries and get the
            # nodes in the transcript path
            for n in split_exons(t, boundaries):
                ref_node_dict[n][t.strand].append(t)
            # add to introns
            for start,end in t.iterintrons():
                ref_intron_dict[(t.strand, start, end)].append(t)
                all_introns.add((t.strand,start,end))
        else:
            if t.strand != NO_STRAND:
                score = float(t.attrs[GTFAttr.SCORE])
                for n in split_exons(t, boundaries):
                    node_score_dict[n][t.strand] += score
            inp_transcripts.append(t)
            # add to introns
            for start,end in t.iterintrons():
                all_introns.add((t.strand,start,end))
    # index introns for fast intersection
    intron_tree = IntervalTree()
    for strand,start,end in all_introns:
        intron_tree.insert_interval(Interval(start,end,strand=strand))
    del all_introns
    # categorize transcripts
    strand_transcript_lists = [[], [], []]
    for t in inp_transcripts:
        # get transcript nodes and introns
        nodes = list(split_exons(t, boundaries))
        introns = set(t.iterintrons())
        # try to resolve strand
        strand = t.strand
        if strand == NO_STRAND:
            strand = resolve_strand(nodes, node_score_dict, ref_node_dict)
        # define opposite strand
        if strand == NO_STRAND:
            opp_strand = NO_STRAND
        else:
            opp_strand = (strand + 1) % 2
        # get all reference transcripts that share introns
        intron_ref_dict = {}
        for start,end in introns:
            if (strand, start, end) in ref_intron_dict:
                refs = ref_intron_dict[(strand, start, end)]
                intron_ref_dict.update((ref.attrs[GTFAttr.TRANSCRIPT_ID],ref) 
                                       for ref in refs)
        intron_refs = []
        for ref in intron_ref_dict.itervalues():
            intron_refs.append((ref,list(split_exons(ref, boundaries))))
        # get all reference transcripts that share coverage
        same_strand_ref_dict = {}
        opp_strand_ref_dict = {}
        for n in nodes:
            if n in ref_node_dict:
                strand_refs = ref_node_dict[n]
                same_strand_ref_dict.update((ref.attrs[GTFAttr.TRANSCRIPT_ID],ref) 
                                            for ref in strand_refs[strand])
                opp_strand_ref_dict.update((ref.attrs[GTFAttr.TRANSCRIPT_ID],ref) 
                                           for ref in strand_refs[opp_strand])
        same_strand_refs = []
        for ref in same_strand_ref_dict.itervalues():
            same_strand_refs.append((ref,list(split_exons(ref, boundaries))))
        opp_strand_refs = []
        for ref in opp_strand_ref_dict.itervalues():            
            opp_strand_refs.append((ref,list(split_exons(ref, boundaries))))
        # categorize
        cinf = categorize_transcript(t, nodes, introns, 
                                     intron_refs,
                                     same_strand_refs,
                                     opp_strand_refs,
                                     intron_tree,
                                     ignore_test=False)
        if cinf.is_test:
            # recategorize test transcripts
            cinf2 = categorize_transcript(t, nodes, introns, 
                                          intron_refs,
                                          same_strand_refs,
                                          opp_strand_refs,
                                          intron_tree,
                                          ignore_test=True)
            cinf = cinf._replace(category=cinf2.category)
        # add annotation attributes
        best_ref_id = (cinf.ref.attrs[GTFAttr.TRANSCRIPT_ID] 
                       if cinf.ref is not None else 'na')
        t.attrs[GTFAttr.CATEGORY] = cinf.category
        t.attrs[GTFAttr.TEST] = '1' if cinf.is_test else '0'
        t.attrs[GTFAttr.ANN_REF_ID] = best_ref_id
        t.attrs[GTFAttr.ANN_COV_RATIO] = cinf.ann_cov_ratio
        t.attrs[GTFAttr.ANN_INTRON_RATIO] = cinf.ann_intron_ratio
        # group transcripts by strand
        strand_transcript_lists[strand].append(t)
    # explictly delete large data structures
    del ref_intron_dict
    del ref_node_dict
    del node_score_dict
    del intron_tree
    del inp_transcripts
    # annotate score and recurrence for transcripts
    for strand_transcripts in strand_transcript_lists:
        # find the intron domains of the transcripts
        boundaries = find_exon_boundaries(strand_transcripts)
        # gather node score/recurrence data
        new_data_func = lambda: {'ids': set(), 
                                 'score': 0.0, 
                                 'pct': 0.0}
        node_data = collections.defaultdict(new_data_func)
        for t in strand_transcripts:
            sample_id = t.attrs[gtf_sample_attr]
            score = float(t.attrs[GTFAttr.SCORE])
            pctrank = float(t.attrs[GTFAttr.PCTRANK])
            # split exons that cross boundaries and to get the
            # nodes in the transcript path
            for n in split_exons(t, boundaries):
                nd = node_data[n]
                nd['ids'].add(sample_id)
                nd['score'] += score
                nd['pct'] += pctrank
        # calculate recurrence and score statistics
        for t in strand_transcripts:
            nodes = list(split_exons(t, boundaries))
            mean_score, mean_pctrank, mean_recur = \
                compute_recurrence_and_score(nodes, node_data)
            t.attrs[GTFAttr.MEAN_SCORE] = mean_score
            t.attrs[GTFAttr.MEAN_PCTRANK] = mean_pctrank
            t.attrs[GTFAttr.MEAN_RECURRENCE] = mean_recur

def annotate_gtf_worker(input_queue, gtf_file, gtf_sample_attr): 
    fileh = open(gtf_file, 'w')
    while True:
        lines = input_queue.get()
        if len(lines) == 0:
            break             
        transcripts = transcripts_from_gtf_lines(lines)
        annotate_locus(transcripts, gtf_sample_attr) 
        for t in transcripts:
            for f in t.to_gtf_features():
                print >>fileh, str(f)
        input_queue.task_done()
        # explicitly delete large objects
        del lines
        del transcripts
    fileh.close()
    input_queue.task_done()

def annotate_gtf_parallel(input_gtf_file,
                          output_gtf_file, 
                          gtf_sample_attr, 
                          num_processors, 
                          tmp_dir):
    # create queue
    input_queue = JoinableQueue(maxsize=num_processors*3)
    # start worker processes
    procs = []
    worker_gtf_files = []
    for i in xrange(num_processors):
        worker_gtf_file = os.path.join(tmp_dir, "annotate_worker%03d.gtf" % (i))
        worker_gtf_files.append(worker_gtf_file)
        args = (input_queue, worker_gtf_file, gtf_sample_attr)
        p = Process(target=annotate_gtf_worker, args=args)
        p.daemon = True
        p.start()
        procs.append(p)
    for lines in parse_loci(open(input_gtf_file)):
        input_queue.put(lines)
    # stop workers
    for p in procs:
        input_queue.put([])
    # close queue
    input_queue.join()
    input_queue.close()
    # join worker processes
    for p in procs:
        p.join()
    # merge/sort worker gtf files
    logging.debug("Merging %d worker GTF file(s)" % (num_processors))
    merge_sort_gtf_files(worker_gtf_files, output_gtf_file, tmp_dir=tmp_dir)
    # remove worker gtf files
    for filename in worker_gtf_files:
        if os.path.exists(filename):
            os.remove(filename)

def main():
    # parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("-v", "--verbose", action="store_true", 
                        dest="verbose", default=False)
    parser.add_argument("-p", type=int, dest="num_processors", default=1,
                        help="Number of processes to run in parallel "
                        "[default=%(default)s]")
    parser.add_argument("--gtf-sample-attr", dest="gtf_sample_attr", 
                        default=GTFAttr.SAMPLE_ID, metavar="ATTR",
                        help="GTF attribute field used to distinguish "
                        "independent samples in order to compute "
                        "recurrence [default=%(default)s]")
    parser.add_argument("run_dir")
    args = parser.parse_args()
    # check command line parameters
    if not os.path.exists(args.run_dir):
        parser.error("Run directory %s not found" % (args.run_dir))
    num_processors = max(1, args.num_processors)
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
    logging.info("num processors:       %d" % (args.num_processors))
    logging.info("gtf sample attribute: %s" % (args.gtf_sample_attr))
    logging.info("run directory:        %s" % (args.run_dir))
    logging.info("----------------------------------")   
    # setup results
    results = config.AssemblylineResults(args.run_dir)
    # function to gather transcript attributes
    logging.info("Annotating GTF file")
    annotate_gtf_parallel(results.transcripts_gtf_file,
                          results.annotated_transcripts_gtf_file,
                          args.gtf_sample_attr,
                          num_processors,
                          results.tmp_dir)
    logging.info("Done")
    return 0

if __name__ == "__main__":
    sys.exit(main())
