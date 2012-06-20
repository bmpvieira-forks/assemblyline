'''
Created on Dec 7, 2011

@author: mkiyer
'''
import operator
import networkx as nx

from transcript import Exon
from assembler_base import TRANSCRIPT_IDS, NODE_DENSITY, NODE_LENGTH

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
            if 'chain' in d:
                new_nodes.extend(d['chain'])
            if 'chain_data' in d:
                chain_data.update(d['chain_data'])
        new_nodes.extend(nodes)
        # sort nodes by genome position and find the min/max
        sorted_nodes = sorted(new_nodes, key=operator.attrgetter('start'))
        # add node attributes
        chain_data.update(dict((n, G.node[n]) for n in nodes))
        attr_dict = {'chain': sorted_nodes,
                     'chain_data': chain_data,
                     'chain_edges': []} 
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
            H.node[u_chain_node]['chain_edges'].append((u,v,d))
    return H

def can_collapse(G,u,v):
    # see if edge nodes have degree larger than '1'
    if ((G.out_degree(u) > 1) or (G.in_degree(v) > 1)):
        return False
    return True

def get_chains(G):
    """
    group nodes into chains
    
    returns a dict mapping node -> chain, as well as a 
    dict mapping chains to nodes
    """
    imin2 = lambda x,y: x if x<=y else y 
    imax2 = lambda x,y: x if x>=y else y 
    node_chain_map = {}
    chains = {}
    # initialize each node to be in a "chain" by itself
    for n in G.nodes_iter():
        node_chain_map[n] = n
        chains[n] = set((n,))
    for u,v in G.edges_iter():
        if not can_collapse(G,u,v):
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
    return node_chain_map, chains

def recalc_strand_specific_graph_attributes(G):
    """
    computes density, length, and transcript ids after graph
    is divided into strand-specific subgraphs and collapsed
    """
    for n,d in G.nodes_iter(data=True):
        chain_nodes = d['chain']
        chain_data = d['chain_data']
        # calculate density of all nodes in chain
        total_mass = 0.0
        total_length = 0
        ids = set()
        for cn in chain_nodes:
            ids.update(chain_data[cn][TRANSCRIPT_IDS])
            length = (cn.end - cn.start)
            total_mass += length * chain_data[cn][NODE_DENSITY]
            total_length += length
        density = total_mass / float(total_length)
        # set attributes
        d[TRANSCRIPT_IDS] = ids
        d[NODE_DENSITY] = density
        d[NODE_LENGTH] = total_length

def collapse_strand_specific_graph(G):
    """
    find groups of nodes that have a single path through them
    and merges them into chains
    
    NOTE: assumes a strand-specific graph.
  
    returns new DiGraph object.  each node has a 'chain' attribute 
    containing the child nodes making up the chain.  nodes also 
    have 'chain_data' and 'chain_edges' attributes with node 
    attribute data and edge data of child nodes  
    """
    node_chain_map, chains = get_chains(G)
    H = add_chains(G, chains, node_chain_map)
    recalc_strand_specific_graph_attributes(H)
    return H
