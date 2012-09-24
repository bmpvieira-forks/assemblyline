'''
Created on Dec 17, 2011

@author: mkiyer
'''
import logging
import bisect

from assemblyline.lib.transcript import POS_STRAND, NEG_STRAND, NO_STRAND, strand_int_to_str
from base import NODE_SCORE
from collapse import get_chains

def find_intron_starts_and_ends(transcripts):
    '''
    input: list of Transcript objects
    output: intron starts/ends by strand (list to index strand, list of starts)    
    
    list of lists of intron starts, list of lists of intron ends
    '''
    # TODO: combine 'find_exon_boundaries' with this
    intron_starts = [set(), set(), set()]
    intron_ends = [set(), set(), set()]
    for transcript in transcripts:
        exons = transcript.exons
        strand = transcript.strand
        if len(exons) == 1:
            continue
        if transcript.strand == NO_STRAND:
            logging.warning("Multi-exonic unstranded transcript detected!")
            continue
        intron_starts[strand].add(exons[0].end)
        for e in transcript.exons[1:-1]:
            intron_starts[strand].add(e.end)
            intron_ends[strand].add(e.start)
        intron_ends[strand].add(exons[-1].start)
    # combine pos/neg to get unstranded
    intron_starts[NO_STRAND] = intron_starts[POS_STRAND].union(intron_starts[NEG_STRAND])
    intron_ends[NO_STRAND] = intron_ends[POS_STRAND].union(intron_ends[NEG_STRAND])
    # sort lists
    for strand in (POS_STRAND, NEG_STRAND, NO_STRAND):
        intron_starts[strand] = sorted(intron_starts[strand])
        intron_ends[strand] = sorted(intron_ends[strand])
    return intron_starts, intron_ends

def trim_left(start, end, intron_starts, intron_ends, overhang_threshold):
    trim_start = start
    # search for the nearest intron end greater than
    # or equal to the current exon start
    i = bisect.bisect_left(intron_ends, start)
    if i != len(intron_ends):
        nearest_intron_end = intron_ends[i]
        # cannot trim past end of exon and cannot
        # trim more than the overhang threshold
        if ((nearest_intron_end < end) and
            (nearest_intron_end <= start + overhang_threshold)):
            trim_start = nearest_intron_end
    return trim_start

def trim_right(start, end, intron_starts, intron_ends, overhang_threshold):
    trim_end = end
    if trim_right:
        # search for nearest intron start less than
        # or equal to the current exon end
        i = bisect.bisect_right(intron_starts, end)
        if i > 0:
            nearest_intron_start = intron_starts[i-1]
            # cannot trim past start of exon and cannot
            # trim more than the overhang threshold
            if ((nearest_intron_start > start) and
                (nearest_intron_start >= end - overhang_threshold)):
                trim_end = nearest_intron_start    
    return trim_end

def trim_transcript(exons, intron_starts, intron_ends, 
                    overhang_threshold):
    # only the very first and last exon of each transcript can be 
    # trimmed
    trim_start = trim_left(exons[0].start, exons[0].end, 
                           intron_starts, intron_ends, 
                           overhang_threshold)
    trim_end = trim_right(exons[-1].start, exons[-1].end, 
                          intron_starts, intron_ends, 
                          overhang_threshold)
    # if both start and end were trimmed it is
    # possible that they could be trimmed to 
    # match different introns and generate a 
    # weird and useless trimmed exon that is not
    # compatible with either the left or right
    # introns. resolve these situations by choosing
    # the smaller of the two trimmed distances as
    # the more likely
    if trim_end <= trim_start:
        left_trim_dist = trim_start - exons[0].start
        right_trim_dist = exons[-1].end - trim_end
        assert left_trim_dist >= 0
        assert right_trim_dist >= 0
        if left_trim_dist <= right_trim_dist:
            trim_start, trim_end = trim_start, exons[-1].end
        else:
            trim_start, trim_end = exons[0].start, trim_end
    assert trim_start < trim_end
    return trim_start, trim_end

def trim_transcripts(transcripts, overhang_threshold):
    if overhang_threshold <= 0:
        return
    # find the set of intron boundaries that govern 
    # trimming
    intron_starts, intron_ends = find_intron_starts_and_ends(transcripts)    
    for t in transcripts:
        trim_start, trim_end = trim_transcript(t.exons, 
                                               intron_starts[t.strand],
                                               intron_ends[t.strand],
                                               overhang_threshold)
        assert trim_start < trim_end
        # finally, modify the transcript
        t.exons[0].start = trim_start
        t.exons[-1].end = trim_end
        

def _trim_utr_nodes(G, nodes, coverage_fraction):
    """
    return list of nodes that should be clipped
    """
    before_total_score = 0
    after_total_score = sum(G.node[nodes[j]][NODE_SCORE] for j in xrange(len(nodes)))
    trim_index = 0
    for i in xrange(1,len(nodes)):
        score = G.node[nodes[i-1]][NODE_SCORE]        
        before_total_score += score
        after_total_score -= score
        before_avg_score = before_total_score / float(i)
        after_avg_score = after_total_score / float(len(nodes) - i)
        fraction = before_avg_score / after_avg_score
        if fraction <= coverage_fraction:
            trim_index = i
    trim_nodes = []
    if trim_index > 0:
        trim_nodes = nodes[:trim_index] 
    return trim_nodes, nodes[trim_index]

def get_start_end_nodes(G):
    start_nodes = []
    end_nodes = []
    both_nodes = []
    for n in G.nodes_iter():
        in_degree = G.in_degree(n)
        out_degree = G.out_degree(n)
        if (in_degree == 0) and (out_degree == 0):
            both_nodes.append(n)
        elif G.in_degree(n) == 0:
            print '5UTR', n
            start_nodes.append(n)
        elif G.out_degree(n) == 0:
            print '3UTR', n
            end_nodes.append(n)
    return start_nodes, end_nodes, both_nodes

def trim_utrs(G, coverage_fraction):
    start_nodes, end_nodes, both_nodes = get_start_end_nodes(G)
    nodes_to_remove = set()
    for n in start_nodes:
        nodes = [n]
        # extend forward
        succs = G.successors(nodes[-1])
        while (len(succs) == 1) and (G.in_degree(succs[0]) == 1):
            nodes.append(succs[0])
            succs = G.successors(nodes[-1])
        # TODO: single exon transcripts will break
        # because trimming could happen from both ends
        if (len(succs) < 1) or (len(nodes) == 1):
            continue
        trim_nodes, new_utr_node = _trim_utr_nodes(G, nodes, coverage_fraction)
        nodes_to_remove.update(trim_nodes)
    for n in end_nodes:
        nodes = [n]
        # extend backward
        preds = G.predecessors(nodes[-1])
        while (len(preds) == 1) and (G.out_degree(preds[0]) == 1):
            nodes.append(preds[0])
            preds = G.predecessors(nodes[-1])
        # TODO: single exon transcripts will break
        # because trimming could happen from both ends
        if (len(preds) < 1) or (len(nodes) == 1):
            continue
        trim_nodes, new_utr_node = _trim_utr_nodes(G, nodes, coverage_fraction)
        nodes_to_remove.update(trim_nodes)
    # remove nodes from graph
    G.remove_nodes_from(nodes_to_remove)
    logging.debug("\t\tTrimmed %d nodes from UTRs" % (len(nodes_to_remove)))

def get_introns(G, reverse=False):
    '''
    input: DiGraph G
    output: introns as (start,end) tuples
    '''
    if reverse:
        G.reverse(copy=False)
    introns = set()
    for u,nbrdict in G.adjacency_iter():
        for v in nbrdict:
            if u.end == v.start:
                continue
            introns.add((u.end,v.start))
    if reverse:
        G.reverse(copy=False)
    return introns

def get_intron_starts_and_ends(G, reverse=False):
    '''
    input: DiGraph G
    output: introns as (start,end) tuples
    '''
    if reverse:
        G.reverse(copy=False)
    intron_starts = set()
    intron_ends = set()
    for u,nbrdict in G.adjacency_iter():
        for v in nbrdict:
            if u.end == v.start:
                continue
            intron_starts.add(u.end)
            intron_ends.add(v.start)
    if reverse:
        G.reverse(copy=False)
    return intron_starts, intron_ends
    
def trim_intron(G, chains, node_chain_map, 
                predecessor_dict, successor_dict, 
                parent, nodes, coverage_fraction):
    # find the nodes immediately upstream of the intron
    pred = None
    for pred in predecessor_dict[nodes[0]]:
        if (parent.end == pred.start) or (parent.start == pred.end):
            break
    if pred is None:
        logging.debug("%s %s" % (str(nodes[-1]), str(successor_dict[nodes[-1]])))
    assert pred is not None
    # calculate the average score of the nodes bordering the intron
    pred_chain_node = node_chain_map[pred]
    pred_nodes = chains[pred_chain_node]
    pred_total_score = sum(G.node[n][NODE_SCORE] for n in pred_nodes)
    pred_avg_score = pred_total_score / len(pred_nodes)
    # find the nodes immediately downstream of the intron
    succ = None
    for succ in successor_dict[nodes[-1]]:
        if (parent.end == succ.start) or (parent.start == succ.end):
            break
    if succ is None:
        logging.debug("%s %s" % (str(nodes[-1]), str(successor_dict[nodes[-1]])))
    assert succ is not None
    # calculate the average score of the nodes bordering the intron
    succ_chain_node = node_chain_map[succ]
    succ_nodes = chains[succ_chain_node]
    succ_total_score = sum(G.node[n][NODE_SCORE] for n in succ_nodes)
    succ_avg_score = succ_total_score / len(succ_nodes)    
    # remove intron nodes with score less than the bordering exons
    cutoff_score = coverage_fraction * max(pred_avg_score, succ_avg_score)
    trim_nodes = set()
    for n in nodes:
        score = G.node[n][NODE_SCORE]
        if score < cutoff_score:
            trim_nodes.add(n)
    return trim_nodes
    
def trim_utr(G, nodes, min_trim_length, coverage_fraction):
    """
    return list of nodes that should be clipped
    """
    # establish seed nodes at least 'min_trim_length' long
    seed_end = 1
    seed_score = G.node[nodes[0]][NODE_SCORE]
    seed_length = (nodes[0].end - nodes[0].start)
    while (seed_length < min_trim_length):
        seed_score += G.node[nodes[seed_end]][NODE_SCORE]
        seed_length += (nodes[seed_end].end - nodes[seed_end].start)
        seed_end += 1
    # extend seed until ratio between seed score and 
    # utr score dips below coverage fraction
    trim_score = sum(G.node[nodes[j]][NODE_SCORE] for j in xrange(seed_end, len(nodes)))
    i = seed_end
    while i < len(nodes):
        seed_avg_score = seed_score / float(i)
        trim_avg_score = trim_score / float(len(nodes) - i)
        frac = trim_avg_score / seed_avg_score
        if frac < coverage_fraction:
            break
        score = G.node[nodes[i]][NODE_SCORE]        
        seed_score += score
        trim_score -= score
        i += 1
    return nodes[i:]

def trim_bidirectional(G, nodes, min_trim_length, coverage_fraction):
    # find max node and use as seed
    seed_index = None
    seed_score = None
    for i,n in enumerate(nodes):
        score = G.node[n][NODE_SCORE]
        if (seed_index is None) or (score > seed_score):
            seed_index = i
            seed_score = score
    #print "nodes", nodes
    #print "seed", seed_index, nodes[seed_index], seed_score
    # extend seed nodes until length greater than min_length
    seed_start = seed_index
    seed_end = seed_index
    seed_length = (nodes[seed_index].end - nodes[seed_index].start)
    while ((seed_length < min_trim_length) and
           ((seed_start > 0) or (seed_end < len(nodes)-1))):
        if seed_start == 0:
            pred_score = 0.0
        else:
            pred_score = G.node[nodes[seed_start-1]][NODE_SCORE]
        if seed_end == (len(nodes)-1):
            succ_score = 0.0
        else:
            succ_score = G.node[nodes[seed_end+1]][NODE_SCORE]
        if (succ_score > pred_score):
            seed_end += 1
            seed_length += nodes[seed_end].end - nodes[seed_end].start
            seed_score += succ_score
        else:
            seed_start -= 1
            seed_length += nodes[seed_start].end - nodes[seed_start].start
            seed_score += pred_score
        #print "extend", seed_start, seed_end, seed_length
    # compute seed score and trimming score cutoff
    seed_avg_score = seed_score / float(seed_end - seed_start + 1)
    score_cutoff = coverage_fraction * seed_avg_score
    trim_nodes = []
    # trim left
    if seed_start > 0:
        for i in xrange(seed_start-1, -1, -1):
            #print nodes[i], G.node[nodes[i]][NODE_SCORE], score_cutoff
            if G.node[nodes[i]][NODE_SCORE] < score_cutoff:
                trim_nodes.extend(nodes[:i+1])
                break
    # trim_right
    if (seed_end+1) < len(nodes):
        for i in xrange(seed_end+1, len(nodes)):
            #print nodes[i], G.node[nodes[i]][NODE_SCORE], score_cutoff
            if G.node[nodes[i]][NODE_SCORE] < score_cutoff:
                trim_nodes.extend(nodes[i:])
                break
    #print "trim", trim_nodes
    return trim_nodes

def trim_graph(G, strand,
               min_trim_length, 
               trim_utr_fraction,
               trim_intron_fraction):
    # get introns
    introns = get_introns(G, reverse=(strand == NEG_STRAND))
    # get lookup table of node predecessors and successors
    successor_dict = {}
    for n,nbrdict in G.adjacency_iter():
        successor_dict[n] = nbrdict.keys()
    predecessor_dict = {}
    G.reverse(copy=False)
    for n,nbrdict in G.adjacency_iter():
        predecessor_dict[n] = nbrdict.keys()
    G.reverse(copy=False)        
    # get 'chains' of contiguous nodes with edge degree of one or less
    node_chain_map, chains = get_chains(G)
    # trim chains
    all_trim_nodes = set()
    for parent, nodes in chains.iteritems():
        trim_nodes = []
        if strand == NEG_STRAND:
            nodes.reverse()
        in_degree = len(predecessor_dict[nodes[0]])
        out_degree = len(successor_dict[nodes[-1]])
        if ((in_degree == 1) and (out_degree == 1) and
            (parent.start, parent.end) in introns): 
            trim_nodes = trim_intron(G, chains, node_chain_map, 
                                     predecessor_dict, successor_dict, 
                                     parent, nodes, trim_intron_fraction)
        else:
            if (in_degree == 0) and (out_degree == 0):
                trim_nodes = trim_bidirectional(G, nodes, min_trim_length, trim_utr_fraction)
            elif in_degree == 0:
                trim_nodes = trim_utr(G, nodes[::-1], min_trim_length, trim_utr_fraction)
            elif out_degree == 0:
                trim_nodes = trim_utr(G, nodes, min_trim_length, trim_utr_fraction)
        all_trim_nodes.update(trim_nodes)
    G.remove_nodes_from(all_trim_nodes)
    if len(all_trim_nodes) > 0:
        logging.debug("\t\t(%s) trimmed %d nodes from graph" % (strand_int_to_str(strand), len(all_trim_nodes)))
