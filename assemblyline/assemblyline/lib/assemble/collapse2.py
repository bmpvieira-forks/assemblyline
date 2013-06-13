'''
Created on Dec 7, 2011

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
import operator
import networkx as nx

from base import Exon, TRANSCRIPT_IDS, NODE_SCORE, NODE_LENGTH, \
    CHAIN_NODES, CHAIN_EDGES, CHAIN_DATA

def can_collapse(G,u,v):
    # see if edge nodes have degree larger than '1'
    if ((G.out_degree(u) > 1) or (G.in_degree(v) > 1)):
        return False
    return True

def can_collapse_contiguous(G,u,v):
    '''
    same as can_collapse but does not allow introns to be collapsed
    '''
    # see if edge nodes have degree larger than '1'
    if ((G.out_degree(u) > 1) or (G.in_degree(v) > 1)):
        return False
    # collapse non-intron edges    
    if (u.end == v.start) or (v.end == u.start):
        return True
    return False

def get_chains(G, introns=True):
    """
    group nodes into chains
    
    returns a dict mapping node -> chain, as well as a 
    dict mapping chains to nodes
    """
    if introns:
        can_collapse_func = can_collapse
    else:
        can_collapse_func = can_collapse_contiguous
    imin2 = lambda x,y: x if x<=y else y 
    imax2 = lambda x,y: x if x>=y else y 
    node_chain_map = {}
    chains = {}
    # initialize each node to be in a "chain" by itself
    for n in G.nodes_iter():
        node_chain_map[n] = n
        chains[n] = set((n,))
    for u,v in G.edges_iter():
        if not can_collapse_func(G,u,v):
            continue
        # get chains containing these nodes
        u_new = node_chain_map[u]
        u_chain = chains[u_new]
        del chains[u_new]
        v_new = node_chain_map[v]
        v_chain = chains[v_new]
        del chains[v_new]
        # merge chains        
        merged_chain = u_chain.union(v_chain)
        merged_node = Exon(imin2(u_new.start, v_new.start),
                           imax2(u_new.end, v_new.end))
        # point all nodes in chain to new parent
        for n in merged_chain:
            node_chain_map[n] = merged_node
        chains[merged_node] = merged_chain
    # sort chain nodes by genome position and store as list
    for parent in chains:
        chains[parent] = sorted(chains[parent], key=operator.attrgetter('start'))
    return node_chain_map, chains

def add_chains(G, chains, node_chain_map):
    H = nx.DiGraph()
    # add chain nodes
    for parent, nodes in chains.iteritems():
        # if nodes have already been collapsed they
        # might have chain attributes and these must
        # be preserved
        new_nodes = []
        chain_data = {}
        for n in nodes:
            d = G.node[n]
            if CHAIN_NODES in d:
                new_nodes.extend(d[CHAIN_NODES])
            if CHAIN_DATA in d:
                chain_data.update(d[CHAIN_DATA])
        new_nodes.extend(nodes)
        # sort nodes by genome position and find the min/max
        sorted_nodes = sorted(new_nodes, key=operator.attrgetter('start'))
        # add node attributes
        chain_data.update(dict((n, G.node[n]) for n in nodes))
        attr_dict = {CHAIN_NODES: sorted_nodes,
                     CHAIN_DATA: chain_data,
                     CHAIN_EDGES: []} 
        H.add_node(parent, attr_dict=attr_dict)
    # add chain edges
    for u,v,d in G.edges_iter(data=True):
        u_chain_node = node_chain_map[u]
        v_chain_node = node_chain_map[v]
        if u_chain_node != v_chain_node:
            H.add_edge(u_chain_node, v_chain_node, attr_dict=d)
        else:
            # add internal chain edge attribute data
            # as 'chain_edges' attribute of parent node
            H.node[u_chain_node][CHAIN_EDGES].append((u,v,d))
    return H

def recalc_strand_specific_graph_attributes(G):
    """
    computes score, length, and transcript ids after graph
    is divided into strand-specific subgraphs and collapsed
    """
    for n,d in G.nodes_iter(data=True):
        chain_nodes = d[CHAIN_NODES]
        chain_data = d[CHAIN_DATA]
        total_length = 0
        max_score = 0.0
        for cn in chain_nodes:
            total_length += (cn.end - cn.start)
            score = chain_data[cn][NODE_SCORE]
            if score > max_score:
                max_score = score
        # set attributes
        d[NODE_LENGTH] = total_length       
        d[NODE_SCORE] = max_score

def collapse_strand_specific_graph(G, introns=True):
    """
    find groups of nodes that have a single path through them
    and merges them into chains
    
    NOTE: assumes a strand-specific graph.
  
    returns new DiGraph object.  each node has a 'chain' attribute 
    containing the child nodes making up the chain.  nodes also 
    have 'chain_data' and 'chain_edges' attributes with node 
    attribute data and edge data of child nodes  
    """
    # TODO: may not need transcript map here
    node_chain_map, chains = get_chains(G, introns)
    H = add_chains(G, chains, node_chain_map)
    recalc_strand_specific_graph_attributes(H)
    return H, node_chain_map
