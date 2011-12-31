'''
Created on Dec 9, 2011

@author: mkiyer
'''
import logging
import collections
import operator
import networkx as nx
import numpy as np

from assembler_base import NODE_DENSITY, NODE_LENGTH, EDGE_OUT_FRAC, EDGE_IN_FRAC, \
    MIN_NODE_DENSITY
from base import imax2, imin2

# for dynamic programming algorithm
TMP_NODE_DENSITY = 'tmpnd'
TMP_EDGE_OUT_FRAC = 'tmpoutfrac'
TMP_EDGE_IN_FRAC = 'tmpinfrac'

PATH_NUM_NODES = 'pn'
PATH_DENSITY_INVSUM = 'pdinvsum'
PATH_DENSITY_MIN = 'pdmin'
PATH_PREV = 'pprev'

def calculate_edge_fractions(G, density_attr=NODE_DENSITY,
                             in_frac_attr=EDGE_IN_FRAC,
                             out_frac_attr=EDGE_OUT_FRAC):
    for u,nbrdict in G.adjacency_iter():
        # find total coverage flowing out of this node
        out_density_total = sum(G.node[v][density_attr] for v in nbrdict)
        for v, eattrs in nbrdict.iteritems():
            # find fraction flowing out of each node
            out_frac = G.node[v][density_attr]/float(out_density_total)
            eattrs[out_frac_attr] = out_frac            
    # now fraction of weight flowing in to each node
    # reverse edge directions and use adjacency iter again
    G.reverse(copy=False)
    for u,nbrdict in G.adjacency_iter():
        # find total coverage density flowing into this node
        in_density_total = sum(G.node[v][density_attr] for v in nbrdict)
        for v, eattrs in nbrdict.iteritems():
            # find fraction flowing into node            
            in_frac = G.node[v][density_attr]/float(in_density_total)
            eattrs[in_frac_attr] = in_frac
    # reverse the edges back to normal
    G.reverse(copy=False)

def init_tmp_attributes(G):
    '''
    set dynamic programming node attributes that are 
    added to the graph temporarily
    '''
    # copy the node weights to a temporary variable so that we can
    # manipulation them in the algorithm and create path attributes
    # for dynamic programming
    for n,d in G.nodes_iter(data=True):
        d[TMP_NODE_DENSITY] = d[NODE_DENSITY]
        d[PATH_NUM_NODES] = 0
        d[PATH_DENSITY_INVSUM] = 0.0
        d[PATH_DENSITY_MIN] = 0.0
        d[PATH_PREV] = None
    for u,nbrdict in G.adjacency_iter():
        for v, d in nbrdict.iteritems():
            d[TMP_EDGE_OUT_FRAC] = d[EDGE_OUT_FRAC]
            d[TMP_EDGE_IN_FRAC] = d[EDGE_IN_FRAC]

def clear_tmp_attributes(G):
    '''
    remove dynamic programming node attributes that are 
    added to the graph temporarily
    '''
    for n,d in G.nodes_iter(data=True):
        del d[TMP_NODE_DENSITY]
        del d[PATH_NUM_NODES]
        del d[PATH_DENSITY_INVSUM]
        del d[PATH_DENSITY_MIN]
        del d[PATH_PREV]
    for u,nbrdict in G.adjacency_iter():
        for v, d in nbrdict.iteritems():
            del d[TMP_EDGE_OUT_FRAC]
            del d[TMP_EDGE_IN_FRAC]

def reset_path_attributes(G):
    """
    must call this before calling the dynamic programming
    algorithm
    """
    # reset path attributes
    for n,d in G.nodes_iter(data=True):
        d[PATH_NUM_NODES] = 0
        d[PATH_DENSITY_INVSUM] = 0.0
        d[PATH_DENSITY_MIN] = 0.0
        d[PATH_PREV] = None

def dynprog_search(G, source, sink):
    """
    Find the highest scoring path by dynamic programming    
    # Taken from NetworkX source code http://networkx.lanl.gov    
    """
    # setup initial path attributes
    reset_path_attributes(G)
    # topological sort allows each node to be visited exactly once
    queue = collections.deque(nx.topological_sort(G))
    while queue:
        u = queue.popleft()
        u_attrs = G.node[u]
        for v in G.successors_iter(u):
            v_attrs = G.node[v]
            # compute path attributes
            if v == sink:
                new_density_invsum = u_attrs[PATH_DENSITY_INVSUM]
                new_num_nodes = u_attrs[PATH_NUM_NODES]
                new_min = u_attrs[PATH_DENSITY_MIN]
            else:
                new_density_invsum = u_attrs[PATH_DENSITY_INVSUM] + (1./v_attrs[TMP_NODE_DENSITY])
                new_num_nodes = u_attrs[PATH_NUM_NODES] + 1
                new_min = imin2(u_attrs[PATH_DENSITY_MIN], v_attrs[TMP_NODE_DENSITY])
            new_hmean = new_num_nodes / new_density_invsum
            # get existing attributes
            if v_attrs[PATH_PREV] is None:
                old_min = 0.0 
                old_hmean = 0.0
            else:
                old_min = v_attrs[PATH_DENSITY_MIN]
                old_hmean = v_attrs[PATH_NUM_NODES] / v_attrs[PATH_DENSITY_INVSUM]
            # update if necessary
            if (new_min, new_hmean) > (old_min, old_hmean):
                v_attrs[PATH_NUM_NODES] = new_num_nodes
                v_attrs[PATH_DENSITY_MIN] = new_min
                v_attrs[PATH_DENSITY_INVSUM] = new_density_invsum
                v_attrs[PATH_PREV] = u

def traceback(G, sink):
    """
    compute path and its density
    """
    path = [sink]
    prev = G.node[sink][PATH_PREV]
    while prev is not None:
        path.append(prev)
        prev = G.node[prev][PATH_PREV]
    path.reverse()
    # find minimum density node along path accounting for edge fractions
    min_density = None
    for i in xrange(1, len(path)-1):        
        u = path[i-1]
        v = path[i]
        w = path[i+1]        
        d = G.node[v]
        # compute length and weight up to this point along path
        v_density_used = (G[u][v][TMP_EDGE_IN_FRAC] * 
                          G[v][w][TMP_EDGE_OUT_FRAC] * 
                          d[TMP_NODE_DENSITY])
        if (min_density is None) or (v_density_used < min_density):
            min_density = v_density_used
    return tuple(path), min_density
                
def find_path(G, source, sink):
    """
    G - graph
    source, sink - start/end nodes
    """    
    # dynamic programming search for best path
    dynprog_search(G, source, sink)
    # traceback to get path
    path, density = traceback(G, sink)
    return path, density

def subtract_path(G, path, density):
    """
    subtract density from nodes along path 
    """
    # first and last node are 'dummy' nodes so don't need to
    # include them in the loop
    for i in xrange(1, len(path)-1):
        u = path[i]
        d = G.node[u]
        d[TMP_NODE_DENSITY] = imax2(MIN_NODE_DENSITY, 
                                    d[TMP_NODE_DENSITY] - density)

def get_adj_matrix_info(G):
    # this works but is numerically unstable
    #A = (np.identity(len(G)) - nx.adj_matrix(G)).I
    A = nx.adj_matrix(G) + np.identity(len(G))
    Y = nx.adj_matrix(G) + np.identity(len(G))
    Y = Y*Y
    while np.sum(Y) > 0:
        A += Y
        A[A > 0] = 1
        Y = Y*Y
    nodes = G.nodes()
    node_indexes = dict((n,i) for i,n in enumerate(nodes))
    return A, nodes, node_indexes

def get_seed_subgraph(H, A, nodes, node_indexes, seed):
    seed_index = node_indexes[seed]
    subgraph_indexes = (A[:,seed_index].getA1() + 
                        A[seed_index,:].getA1()).nonzero()[0]
    subgraph_nodes = tuple(nodes[i] for i in subgraph_indexes)
    G = H.subgraph(subgraph_nodes)
    return G

def find_suboptimal_paths(G, source, sink, 
                          fraction_major_path=1e-3,
                          max_paths=100):
    """
    finds suboptimal paths through graph G using greedy algorithm
    that finds the highest density path, subtracts the path density 
    from the graph, and repeats
    
    algorithm stops when path density falls below a fraction of the
    highest density path found in the graph, or when 'max_paths'
    iterations have completed
    """
    # setup temporary graph attributes
    init_tmp_attributes(G)
    # get adjacency matrix and ordered node list
    A, nodes, node_index_dict = get_adj_matrix_info(G)
    # store paths in a dictionary in order to avoid redundant paths
    # that arise due to issues in the data
    path_results = collections.OrderedDict()
    # maintain list of nodes sorted by density 
    seed_nodes = G.nodes()
    seed_nodes.remove(source)
    seed_nodes.remove(sink)
    seed_nodes.sort(key=lambda n: G.node[n][TMP_NODE_DENSITY])
    # iterate to find suboptimal paths
    logging.debug("\tFinding suboptimal paths")
    iterations = 1
    highest_density = 0.0
    while (len(seed_nodes) > 0) and (iterations < max_paths):
        # get subgraph implied by highest density seed node and re-sort
        seed = seed_nodes.pop()
        Gsub = get_seed_subgraph(G, A, nodes, node_index_dict, seed)
        seed_nodes.sort(key=lambda n: G.node[n][TMP_NODE_DENSITY])
        # find best path and subtract
        path, density = find_path(Gsub, source, sink)
        highest_density = imax2(highest_density, density)
        logging.debug("\t\tdensity=%f" % (density))
        # store path
        if path not in path_results:
            path_results[path] = density
        # subtract path and recalculate edge fractions
        subtract_path(G, path, density)
        calculate_edge_fractions(G, density_attr=TMP_NODE_DENSITY,
                                 in_frac_attr=TMP_EDGE_IN_FRAC,
                                 out_frac_attr=TMP_EDGE_OUT_FRAC)
        iterations +=1
    logging.debug("\t\tpath finding iterations=%d" % iterations)
    # cleanup graph attributes
    clear_tmp_attributes(G) 
    # return (path,density) tuples sorted from high -> low density
    lowest_density = highest_density * fraction_major_path
    path_list = [(p,d) for (p,d) in path_results.iteritems() if d >= lowest_density]
    path_list.sort(key=operator.itemgetter(1), reverse=True)
    return path_list  
    #return sorted(path_results.items(), key=operator.itemgetter(1), reverse=True)

