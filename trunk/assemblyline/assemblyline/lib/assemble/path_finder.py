'''
Created on Dec 9, 2011

@author: mkiyer
'''
import logging
import collections
import networkx as nx

from assemblyline.lib.base import imax2, imin2
from base import NODE_DENSITY, MIN_DENSITY

# for dynamic programming algorithm
TMP_NODE_DENSITY = 'tmpnd'
PATH_MIN_DENSITY = 'pmin'
PATH_PREV = 'pprev'

def init_tmp_attributes(G):
    '''
    set node attributes that are added to the graph temporarily
    '''
    # copy the node weights to a temporary variable so that we can
    # manipulation them in the algorithm and create path attributes
    # for dynamic programming
    for n,d in G.nodes_iter(data=True):
        d[TMP_NODE_DENSITY] = d[NODE_DENSITY]
        d[PATH_MIN_DENSITY] = None
        d[PATH_PREV] = None

def clear_tmp_attributes(G):
    '''
    remove node attributes that are added to the graph temporarily
    '''
    for n,d in G.nodes_iter(data=True):
        del d[TMP_NODE_DENSITY]
        del d[PATH_MIN_DENSITY]
        del d[PATH_PREV]

def reset_path_attributes(G):
    """
    must call this before calling the dynamic programming algorithm
    """
    # reset path attributes
    for n,d in G.nodes_iter(data=True):
        d[PATH_MIN_DENSITY] = None
        d[PATH_PREV] = None

def dynprog_search(G, source):
    """
    Find the highest scoring path by dynamic programming    
    # Adapted from NetworkX source code http://networkx.lanl.gov    
    """
    # setup initial path attributes
    reset_path_attributes(G)
    G.node[source][PATH_MIN_DENSITY] = G.node[source][TMP_NODE_DENSITY] 
    # topological sort allows each node to be visited exactly once
    queue = collections.deque(nx.topological_sort(G))
    while queue:
        u = queue.popleft()
        path_min_density = G.node[u][PATH_MIN_DENSITY]
        for v in G.successors_iter(u):
            v_attrs = G.node[v]
            v_density = v_attrs[TMP_NODE_DENSITY]
            # compute minimum density that would occur if path
            # traversed through node 'v'
            new_min_density = imin2(path_min_density, v_density)
            # update if density is larger
            if ((v_attrs[PATH_PREV] is None) or
                (new_min_density > v_attrs[PATH_MIN_DENSITY])):
                v_attrs[PATH_MIN_DENSITY] = new_min_density
                v_attrs[PATH_PREV] = u

def traceback(G, sink):
    """
    compute path and its density
    """
    path = [sink]
    density = G.node[sink][PATH_MIN_DENSITY]
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
    dynprog_search(G, source)
    # traceback to get path
    path, density = traceback(G, sink)
    return path, density

def subtract_path(G, path, density):
    """
    subtract density from nodes along path 
    """
    for i in xrange(len(path)):
        u = path[i]
        d = G.node[u]
        d[TMP_NODE_DENSITY] = imax2(MIN_DENSITY, 
                                    d[TMP_NODE_DENSITY] - density)

def find_suboptimal_paths(G, fraction_major_path=1e-3, max_paths=1000):
    """
    finds suboptimal paths through graph G using a greedy algorithm that 
    finds the highest density path using dynamic programming, subtracts 
    the path density from the graph, and repeats
    
    paths with density lower than 'fraction_major_path' are not
    returned.  algorithm may stop prematurely if 'max_paths' iterations
    have completed.
    """
    # setup temporary graph attributes
    init_tmp_attributes(G)
    # find source and sink nodes in graph
    sources = set()
    sinks = set()
    for n in G.nodes_iter():
        if G.in_degree(n) == 0:
            sources.add(n)
        if G.out_degree(n) == 0:
            sinks.add(n)
    # TODO: remove assert
    assert len(sources) == 1
    assert len(sinks) == 1
    source = sources.pop()
    sink = sinks.pop()
    # store paths in a dictionary in order to avoid redundant paths
    # that arise when the heuristic assumptions of the algorithm fail
    path_results = collections.OrderedDict()
    # find highest density path
    path, density = find_path(G, source, sink)
    path_results[path] = density
    subtract_path(G, path, density)
    # iterate to find suboptimal paths
    iterations = 1
    highest_density = density
    lowest_density = highest_density * fraction_major_path
    while iterations < max_paths:
        # find path
        path, density = find_path(G, source, sink)
        if density < lowest_density:
            break
        # store path
        if path not in path_results:
            path_results[path] = density
        # TODO: remove assert
        assert highest_density >= density
        # subtract path density from graph and resort seed nodes
        subtract_path(G, path, density)
        iterations +=1
    logging.debug("\t\tpath finding iterations=%d" % iterations)
    # cleanup graph attributes
    clear_tmp_attributes(G) 
    # return (path,density) tuples sorted from high -> low density
    return path_results.items()
