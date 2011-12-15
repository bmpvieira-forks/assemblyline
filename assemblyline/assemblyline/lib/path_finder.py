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
TMP_NODE_DENSITY = 'tmpnw'
TMP_EDGE_OUT_FRAC = 'tmpoutfrac'
TMP_EDGE_IN_FRAC = 'tmpinfrac'

PATH_DENSITY_MINMAX = 'pmd'
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

def init_tmp_attrs(G):
    # copy the node weights to a temporary variable so that we can
    # manipulation them in the algorithm
    for n,d in G.nodes_iter(data=True):
        d[TMP_NODE_DENSITY] = d[NODE_DENSITY]

def clear_tmp_attrs(G):
    for n,d in G.nodes_iter(data=True):
        del d[TMP_NODE_DENSITY]

def init_path_attributes(G):
    """
    must call this before calling the dynamic programming
    algorithm
    """
    # reset path attributes
    for n,d in G.nodes_iter(data=True):
        d[PATH_DENSITY_MINMAX] = (None,None)
        d[PATH_PREV] = None

def clear_path_attributes(G):
    '''
    remove dynamic programming node attributes that are 
    added to the graph temporarily
    
    must call this before calling the dynamic programming
    algorithm again on the same graph    
    '''
    # clear path attributes
    for n,d in G.nodes_iter(data=True):
        del d[PATH_DENSITY_MINMAX]
        del d[PATH_PREV]

def dynprog_search(G, source, sink):
    """
    Find the highest scoring path by dynamic programming    
    # Taken from NetworkX source code http://networkx.lanl.gov    
    """
    def imin2none(x,y):
        if x is None:
            return y
        if y is None:
            return x
        return x if x <= y else y
    def imax2none(x,y):
        if x is None:
            return y
        if y is None:
            return x
        return x if x >= y else y
    # setup initial path attributes
    init_path_attributes(G)
    # topological sort allows each node to be visited exactly once
    queue = collections.deque(nx.topological_sort(G))
    while queue:
        u = queue.popleft()
        minmax = G.node[u][PATH_DENSITY_MINMAX]
        for v in G.successors_iter(u):
            v_attrs = G.node[v]
            if v == sink:
                v_density = None
            else:
                v_density = v_attrs[TMP_NODE_DENSITY]
            # compute new density at node v
            new_minmax = (imin2none(minmax[0], v_density),
                          imax2none(minmax[1], v_density)) 
            # update if density is larger
            if ((v_attrs[PATH_PREV] is None) or
                (new_minmax > v_attrs[PATH_DENSITY_MINMAX])):
                v_attrs[PATH_DENSITY_MINMAX] = new_minmax
                v_attrs[PATH_PREV] = u

def traceback(G, sink):
    """
    compute path and its density
    """
    path = [sink]
    density = G.node[sink][PATH_DENSITY_MINMAX][0]
    prev = G.node[sink][PATH_PREV]
    while prev is not None:
        path.append(prev)
        prev = G.node[prev][PATH_PREV]
    path.reverse()
    return tuple(path), density
                
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
    init_tmp_attrs(G)
    # get adjacency matrix and ordered node list
    A = (np.identity(len(G)) - nx.adj_matrix(G)).I
    nodes = G.nodes()
    node_indexes = dict((n,i) for i,n in enumerate(nodes))
    # store paths in a dictionary in order to avoid redundant paths
    # that arise due to issues in the data
    path_results = collections.OrderedDict()
    # maintain list of nodes sorted by density 
    seed_nodes = G.nodes()
    seed_nodes.sort(key=lambda n: G.node[n][TMP_NODE_DENSITY])
    # get subgraph implied by highest density seed node
    Gsub = get_seed_subgraph(G, A, nodes, node_indexes, seed_nodes[-1])
    # find best path and calculate density lower bound
    logging.debug("\tFinding best path in graph with %d nodes" % len(G)) 
    path, density = find_path(Gsub, source, sink)
    # store path
    if path not in path_results:
        path_results[path] = density
    lowest_density = imax2(MIN_NODE_DENSITY, density * fraction_major_path)
    logging.debug("\t\tdensity=%f lower_bound=%f" % (density, lowest_density))
    # subtract path and resort seed nodes
    subtract_path(G, path, density)
    seed_nodes.sort(key=lambda n: G.node[n][TMP_NODE_DENSITY])
    # iterate to find suboptimal paths
    logging.debug("\tFinding suboptimal paths")
    iterations = 1
    while (iterations < max_paths):
        # get subgraph implied by highest density seed node
        Gsub = get_seed_subgraph(G, A, nodes, node_indexes, seed_nodes[-1])
        # find best path and subtract
        path, density = find_path(Gsub, source, sink)
        if density <= lowest_density:
            break
        logging.debug("\t\tdensity=%f" % (density))
        # store path
        if path not in path_results:
            path_results[path] = density
        # subtract path and resort seed node list
        subtract_path(G, path, density)
        seed_nodes.sort(key=lambda n: G.node[n][TMP_NODE_DENSITY])
        iterations +=1
    logging.debug("\t\tpath finding iterations=%d" % iterations)
    # cleanup graph attributes
    clear_path_attributes(G) 
    clear_tmp_attrs(G)
    # return (path,density) tuples sorted from high -> low density
    return sorted(path_results.items(), key=operator.itemgetter(1), reverse=True)
