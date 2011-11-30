'''
Created on Mar 12, 2011

@author: mkiyer
'''
import collections
import operator
import numpy as np
import networkx as nx

from base import Exon, POS_STRAND, NEG_STRAND, NO_STRAND, NODE_WEIGHT, \
    EDGE_OUT_FRAC

def get_chains(G):
    """find chains of nodes with degree of '1' and smoothen their 
    coverage"""
    node_chain_map = {}
    chains = {}
    last_chain_id = 0
    for n in G.nodes_iter():
        node_chain_map[n] = last_chain_id
        chains[last_chain_id] = set((n,))
        last_chain_id += 1
    for u,v in G.edges_iter():
        # see if edge nodes have degree '1'
        if ((G.out_degree(u) > 1) or (G.in_degree(v) > 1)):
            continue
        # get chains containing these nodes
        u_id = node_chain_map[u]
        u_chain = chains[u_id]
        del chains[u_id]
        v_id = node_chain_map[v]
        v_chain = chains[v_id]
        del chains[v_id]
        # merge chains
        merged_chain = u_chain.union(v_chain)
        merged_id = last_chain_id
        last_chain_id += 1
        node_chain_map[u] = merged_id
        node_chain_map[v] = merged_id
        chains[merged_id] = merged_chain
    return node_chain_map, chains

def add_chain_nodes(H, G, chains):
    chain_id = 0
    for nodes in chains.itervalues():
        total_weight = 0.0
        total_length = 0.0
        tdata = {}
        # gather data from all nodes in chain
        for n in nodes:
            nd = G.node[n]
            total_weight += nd[NODE_WEIGHT]
            total_length += (n.end - n.start)
            if 'data' in nd:
                tdata[n] = nd['data']
        # sort nodes by genome position and find the min/max
        sorted_nodes = sorted(nodes, key=operator.attrgetter('start'))
        # make new node
        H.add_node(chain_id, chain=sorted_nodes, data=tdata, 
                   weight=total_weight, length=total_length)
        chain_id += 1

def add_chain_edges(H, G, node_chain_map):
    for u,v,d in G.edges_iter(data=True):
        u_chain_id = node_chain_map[u]
        v_chain_id = node_chain_map[v]
        if u_chain_id != v_chain_id:
            H.add_edge(u_chain_id, v_chain_id, attr_dict=d)

def collapse_chains(G, strand):
    node_chain_map, chains = get_chains(G)
    H = nx.DiGraph()
    add_chain_nodes(H, G, chains)
    add_chain_edges(H, G, node_chain_map)
    return H

def smoothen_graph_leastsq(G):
    # cannot smooth a graph with just one node
    if len(G) == 1:
        return
    # label nodes
    cur_id = 0
    eqns = {}
    cols = {}
    y = np.zeros(len(G), dtype=np.float)
    # assign numeric ids to nodes
    for n,d in G.nodes_iter(data=True):
        d['id'] = cur_id
        y[cur_id] = d[NODE_WEIGHT]
        if G.in_degree(n) > 0:
            a = np.zeros(len(G), dtype=np.float)
            a[cur_id] = -1.0
            eqns[cur_id] = a
        cols[cur_id] = n
        cur_id += 1
    # setup equations for each node
    for u,v,d in G.edges_iter(data=True):
        u_id = G.node[u]['id']
        v_id = G.node[v]['id']
        out_frac = d[EDGE_OUT_FRAC]
        eqns[v_id][u_id] = out_frac
    # build matrix for least squares error
    # minimization
    row_ids = sorted(eqns)
    A = np.vstack(eqns[id] for id in row_ids)
    A = np.asmatrix(A)
    I = np.identity(len(G), dtype=np.float)
    # perform least squares minimization
    beta = (I - A.T*((A*A.T).I)*A).dot(y)
    # set new node weights
    for col_id,n in cols.iteritems():
        G.node[n][NODE_WEIGHT] = beta[0,col_id]
        del G.node[n]['id']
    #for n,d in G.nodes_iter(data=True):
    #    print 'NEW NODE', n, d



