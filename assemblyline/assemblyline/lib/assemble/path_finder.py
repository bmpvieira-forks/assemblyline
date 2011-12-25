'''
Created on Dec 9, 2011

@author: mkiyer
'''
import logging
import collections
import operator
import networkx as nx
import numpy as np

from base import NODE_DENSITY, MIN_DENSITY, FWD, REV
from assemblyline.lib.base import imax2, imin2

# for dynamic programming algorithm
TMP_NODE_DENSITY = 'tmpnd'

def init_tmp_attributes(G):
    '''
    set node attributes that are 
    added to the graph temporarily
    '''
    # copy the node weights to a temporary variable so that we can
    # manipulation them in the algorithm and create path attributes
    # for dynamic programming
    for n,d in G.nodes_iter(data=True):
        d[TMP_NODE_DENSITY] = d[NODE_DENSITY]

def clear_tmp_attributes(G):
    '''
    remove node attributes that are 
    added to the graph temporarily
    '''
    for n,d in G.nodes_iter(data=True):
        del d[TMP_NODE_DENSITY]

def find_path(G, seed):
    path = collections.deque([seed])
    density = G.node[seed][TMP_NODE_DENSITY]
    while True:
        best_node = None
        best_density = None
        for n in G.successors_iter(path[-1]):
            n_density = G.node[n][TMP_NODE_DENSITY]
            if (best_density is None) or (n_density > best_density):
                best_node = n
                best_density = n_density
        if best_node is None:
            break
        path.append(best_node)
        density = imin2(density, best_density)
    while True:
        best_node = None
        best_density = None
        for n in G.predecessors_iter(path[0]):
            n_density = G.node[n][TMP_NODE_DENSITY]
            if (best_density is None) or (n_density > best_density):
                best_node = n
                best_density = n_density
        if best_node is None:
            break
        path.appendleft(best_node)
        density = imin2(density, best_density)
    return tuple(path), density

def subtract_path(G, path, density):
    """
    subtract density from nodes along path 
    """
    # first and last node are 'dummy' nodes so don't need to
    # include them in the loop
    for i in xrange(1, len(path)-1):
        u = path[i]
        d = G.node[u]
        d[TMP_NODE_DENSITY] = imax2(MIN_DENSITY, d[TMP_NODE_DENSITY] - density)

def find_suboptimal_paths(G, fraction_major_path=1e-3, max_paths=100):
    """
    finds suboptimal paths through graph G using greedy algorithm
    that finds the highest density path, subtracts the path density 
    from the graph, and repeats
    
    paths with density lower than 'fraction_major_path' are not
    returned.  algorithm may stop prematurely if 'max_paths' iterations
    have completed.
    """
    # setup temporary graph attributes
    init_tmp_attributes(G)    
    # store paths in a dictionary in order to avoid redundant paths
    # that arise due to issues in the data
    path_results = collections.OrderedDict()
    # maintain list of nodes sorted by density 
    seed_nodes = G.nodes()
    seed_nodes.sort(key=lambda n: G.node[n][TMP_NODE_DENSITY])    
    # iterate to find suboptimal paths
    logging.debug("\tFinding suboptimal paths in graph with %d nodes" % (len(G)))
    iterations = 0
    highest_density = 0.0
    while (len(seed_nodes) > 0) and (iterations < max_paths):
        # get seed node with highest density and re-sort
        seed = seed_nodes.pop()
        seed_nodes.sort(key=lambda n: G.node[n][TMP_NODE_DENSITY])
        # find path and subtract
        path, density = find_path(G, seed)
        logging.debug("TSS? %s" % ("tss_id" in G.node[path[0]].keys()))
        logging.debug("PATH %s" % (map(str,path)))
        subtract_path(G, path, density)
        highest_density = imax2(highest_density, density)
        logging.debug("\t\tdensity=%f" % (density))
        # store path
        if path not in path_results:
            path_results[path] = density
        iterations +=1
    logging.debug("\t\tpath finding iterations=%d" % iterations)
    # cleanup graph attributes
    clear_tmp_attributes(G) 
    # return (path,density) tuples sorted from high -> low density
    lowest_density = highest_density * fraction_major_path
    path_list = [(p,d) for (p,d) in path_results.iteritems() if d >= lowest_density]
    path_list.sort(key=operator.itemgetter(1), reverse=True)
    return path_list

