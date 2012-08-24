'''
Created on Dec 17, 2011

@author: mkiyer
'''
import logging
import bisect

from assemblyline.lib.transcript import POS_STRAND, NEG_STRAND, NO_STRAND, strand_int_to_str
from base import NODE_DENSITY
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
    before_total_density = 0
    after_total_density = sum(G.node[nodes[j]][NODE_DENSITY] for j in xrange(len(nodes)))
    trim_index = 0
    for i in xrange(1,len(nodes)):
        density = G.node[nodes[i-1]][NODE_DENSITY]        
        before_total_density += density
        after_total_density -= density
        before_avg_density = before_total_density / float(i)
        after_avg_density = after_total_density / float(len(nodes) - i)
        fraction = before_avg_density / after_avg_density
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

def calc_average_density(G, nodes):
    total_mass = 0.0
    total_length = 0
    for n in nodes:
        length = (n.end - n.start)
        total_mass += length * G.node[n][NODE_DENSITY]
        total_length += length
    return total_mass / float(total_length)
    
def trim_intron(G, chains, node_chain_map, successor_dict, 
                parent, nodes, trim_intron_fraction=0.25):
    # find the nodes immediately downstream of the intron
    succ = None
    for succ in successor_dict[nodes[-1]]:
        if (parent.end == succ.start) or (parent.start == succ.end):
            break
    if succ is None:
        logging.debug("%s %s" % (str(nodes[-1]), str(successor_dict[nodes[-1]])))
    assert succ is not None
    # calculate the average density of the nodes bordering the intron
    succ_chain_node = node_chain_map[succ]
    succ_nodes = chains[succ_chain_node]
    succ_density = calc_average_density(G, succ_nodes)
    # remove intron nodes with density less than the bordering exons
    trim_nodes = set()
    for n in nodes:
        density = G.node[n][NODE_DENSITY]
        if density < (trim_intron_fraction * succ_density):
            trim_nodes.add(n)
    return trim_nodes
    
def trim_utr(G, nodes, coverage_fraction):
    """
    return list of nodes that should be clipped
    """
    before_total_density = 0
    after_total_density = sum(G.node[nodes[j]][NODE_DENSITY] for j in xrange(len(nodes)))
    trim_index = 0
    for i in xrange(1,len(nodes)):
        density = G.node[nodes[i-1]][NODE_DENSITY]        
        before_total_density += density
        after_total_density -= density
        before_avg_density = before_total_density / float(i)
        after_avg_density = after_total_density / float(len(nodes) - i)
        fraction = before_avg_density / after_avg_density
        #print i, fraction, before_avg_density, after_avg_density
        if fraction <= coverage_fraction:
            trim_index = i
    trim_nodes = []
    if trim_index > 0:
        trim_nodes = nodes[:trim_index] 
    return trim_nodes
   
def trim_graph(G, strand,
               overhang_threshold, 
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
            trim_nodes = trim_intron(G, chains, node_chain_map, successor_dict, 
                                     parent, nodes, trim_intron_fraction)
        else:
            if (in_degree == 0) and (out_degree == 0):
                # TODO: trim in both directions
                continue
            elif in_degree == 0:
                trim_nodes = trim_utr(G, nodes, trim_utr_fraction)
            elif out_degree == 0:
                trim_nodes = trim_utr(G, nodes[::-1], trim_utr_fraction)
        all_trim_nodes.update(trim_nodes)
    G.remove_nodes_from(all_trim_nodes)
    if len(all_trim_nodes) > 0:
        logging.debug("\t\t(%s) trimmed %d nodes from graph" % (strand_int_to_str(strand), len(all_trim_nodes)))
