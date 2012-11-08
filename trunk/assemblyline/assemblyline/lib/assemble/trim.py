'''
Created on Dec 17, 2011

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

from assemblyline.lib.bx.intersection import Interval, IntervalTree
from assemblyline.lib.transcript import NEG_STRAND, strand_int_to_str

from base import NODE_SCORE
from collapse import get_chains

def get_introns(G, reverse=False):
    '''
    input: DiGraph G
    output: introns as (start,end) tuples
    '''
    if reverse:
        G.reverse(copy=False)
    introns = set()
    intron_tree = IntervalTree()
    for u,nbrdict in G.adjacency_iter():
        for v in nbrdict:
            if u.end == v.start:
                continue
            introns.add((u.end,v.start))
            intron_tree.insert_interval(Interval(u.end,v.start))
    if reverse:
        G.reverse(copy=False)
    return introns, intron_tree

def trim_intron(G, nodes, cutoff_score):
    '''
    remove intron nodes with score less than the bordering exons
    '''
    trim_nodes = set()
    for n in nodes:
        score = G.node[n][NODE_SCORE]
        if score < cutoff_score:
            trim_nodes.add(n)
    return trim_nodes

def trim_intronic_utr(G, nodes, cutoff_score):
    trim_nodes = set()
    for n in nodes:
        if G.node[n][NODE_SCORE] >= cutoff_score:
            break
        trim_nodes.add(n)
    return trim_nodes

#def trim_intron(G, chains, node_chain_map, 
#                predecessor_dict, successor_dict, 
#                parent, nodes, coverage_fraction):
#    # find the nodes immediately upstream of the intron
#    pred = None
#    for pred in predecessor_dict[nodes[0]]:
#        if (parent.end == pred.start) or (parent.start == pred.end):
#            break
#    if pred is None:
#        logging.debug("%s %s" % (str(nodes[-1]), str(successor_dict[nodes[-1]])))
#    assert pred is not None
#    # calculate the average score of the nodes bordering the intron
#    pred_chain_node = node_chain_map[pred]
#    pred_nodes = chains[pred_chain_node]
#    pred_total_score = sum(G.node[n][NODE_SCORE] for n in pred_nodes)
#    pred_avg_score = pred_total_score / len(pred_nodes)
#    # find the nodes immediately downstream of the intron
#    succ = None
#    for succ in successor_dict[nodes[-1]]:
#        if (parent.end == succ.start) or (parent.start == succ.end):
#            break
#    if succ is None:
#        logging.debug("%s %s" % (str(nodes[-1]), str(successor_dict[nodes[-1]])))
#    assert succ is not None
#    # calculate the average score of the nodes bordering the intron
#    succ_chain_node = node_chain_map[succ]
#    succ_nodes = chains[succ_chain_node]
#    succ_total_score = sum(G.node[n][NODE_SCORE] for n in succ_nodes)
#    succ_avg_score = succ_total_score / len(succ_nodes)    
#    # remove intron nodes with score less than the bordering exons
#    cutoff_score = coverage_fraction * max(pred_avg_score, succ_avg_score)
#    trim_nodes = set()
#    for n in nodes:
#        score = G.node[n][NODE_SCORE]
#        if score < cutoff_score:
#            trim_nodes.add(n)
#    return trim_nodes

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

def trim_utr(G, nodes, min_trim_length, coverage_fraction):
    """
    return list of nodes that should be clipped
    """
    # establish seed nodes at least 'min_trim_length' long
    seed_end = 1
    seed_score = G.node[nodes[0]][NODE_SCORE]
    seed_length = (nodes[0].end - nodes[0].start)
    while ((seed_length < min_trim_length) and
           (seed_end < len(nodes))):
        seed_score += G.node[nodes[seed_end]][NODE_SCORE]
        seed_length += (nodes[seed_end].end - nodes[seed_end].start)
        seed_end += 1
    seed_avg_score = seed_score / float(seed_end)
    # find point where trimmed nodes have low score relative to seed
    trim_score = sum(G.node[nodes[j]][NODE_SCORE] for j in xrange(seed_end, len(nodes)))
    i = seed_end
    while i < len(nodes):
        trim_avg_score = trim_score / float(len(nodes) - i)
        frac = trim_avg_score / seed_avg_score
        if frac < coverage_fraction:
            break
        score = G.node[nodes[i]][NODE_SCORE]        
        trim_score -= score
        i += 1
    return nodes[i:]

def trim_graph(G, strand,
               min_trim_length, 
               trim_utr_fraction,
               trim_intron_fraction):
    # get 'chains' of contiguous nodes with edge degree of one or less
    node_chain_map, chains = get_chains(G)
    # setup dictionaries of predecessors and successors
    successor_dict = {}
    for n,nbrdict in G.adjacency_iter():
        successor_dict[n] = nbrdict.keys()
    predecessor_dict = {}
    G.reverse(copy=False)
    for n,nbrdict in G.adjacency_iter():
        predecessor_dict[n] = nbrdict.keys()
    G.reverse(copy=False)
    # setup intron data structures
    introns = {}
    intron_tree = IntervalTree()
    reverse = (strand == NEG_STRAND)
    for u,nbrdict in G.adjacency_iter():
        for v in nbrdict:
            if u.end == v.start:
                continue
            # calculate score of the chains
            u_chain_nodes = chains[node_chain_map[u]]
            u_score = max(G.node[n][NODE_SCORE] for n in u_chain_nodes)
            v_chain_nodes = chains[node_chain_map[v]]
            v_score = max(G.node[n][NODE_SCORE] for n in v_chain_nodes)
            if reverse:
                left, right = v, u
            else:
                left, right = u, v
            introns[(left.end,right.start)] = (u_score, v_score)
            intron_tree.insert_interval(Interval(left.end, right.start, 
                                                 value=(u_score,v_score)))
    # trim chains
    all_trim_nodes = set()
    for parent, nodes in chains.iteritems():
        if strand == NEG_STRAND:
            nodes.reverse()
        in_degree = len(predecessor_dict[nodes[0]])
        out_degree = len(successor_dict[nodes[-1]])
        trim_nodes = set()
        if ((in_degree == 1) and (out_degree == 1) and
            (parent.start, parent.end) in introns): 
            # intron retention - a chain of nodes precisely matches an 
            # intron, so we can potentially remove the entire chain
            pred_score, succ_score = introns[(parent.start, parent.end)]
            cutoff_score = trim_intron_fraction * max(pred_score, succ_score)
            trim_nodes.update(trim_intron(G, nodes, cutoff_score))
        else:
            # determine whether this node chain is intronic. intronic node
            # chains are trimmed more strictly due to intronic pre-mrna
            found_intron = False
            max_pred_score = 0.0
            max_succ_score = 0.0
            for hit in intron_tree.find(parent.start, parent.end):
                # ignore contained introns
                if (hit.start > parent.start) and (hit.end < parent.end):
                    continue
                # set intron flag and keep track of highest coverage 
                # overlapping intron to make trimming conservative
                found_intron = True
                pred_score, succ_score = hit.value
                if pred_score > max_pred_score:
                    max_pred_score = pred_score
                if succ_score > max_succ_score:
                    max_succ_score = succ_score
            if (in_degree == 0) and (out_degree == 0):
                if found_intron:
                    cutoff_score = trim_intron_fraction * max(max_pred_score, max_succ_score)
                    trim_nodes.update(trim_intron(G, nodes, cutoff_score))
                trim_nodes.update(trim_bidirectional(G, nodes, min_trim_length, trim_utr_fraction))
            elif in_degree == 0:
                if found_intron:
                    cutoff_score = trim_intron_fraction * max_succ_score
                    trim_nodes.update(trim_intronic_utr(G, nodes[::-1], cutoff_score))
                trim_nodes.update(trim_utr(G, nodes[::-1], min_trim_length, trim_utr_fraction))
            elif out_degree == 0:
                if found_intron:
                    cutoff_score = trim_intron_fraction * max_pred_score
                    trim_nodes.update(trim_intronic_utr(G, nodes, cutoff_score))
                trim_nodes.update(trim_utr(G, nodes, min_trim_length, trim_utr_fraction))
        all_trim_nodes.update(trim_nodes)
    G.remove_nodes_from(all_trim_nodes)
    if len(all_trim_nodes) > 0:
        logging.debug("\t\t(%s) trimmed %d nodes from graph" % (strand_int_to_str(strand), len(all_trim_nodes)))
