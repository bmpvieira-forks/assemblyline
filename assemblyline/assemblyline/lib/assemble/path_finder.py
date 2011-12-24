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

def get_path_density(G, path, density_attr):
    return min(G.node[n][density_attr] for n in path)

def get_highest_density_kmer(G, kmer_set, density_attr):
    best_density = -1e6
    best_kmer = None
    for kmer in kmer_set:
        density = get_path_density(G, kmer, TMP_NODE_DENSITY)
        if density > best_density:
            best_kmer = kmer
            best_density = density
    return best_kmer, best_density
    
def extend_seed(G, seed_kmers, kmax):
    """
    G - graph
    kmer_dict - path information
    kmax - max length extension to consider
    seed - starting node
    """
    best_kmer = None
    best_density = -1e6
    # TODO: don't use kmax if bordering against a START node
    for k in xrange(kmax,1,-1):        
        kmer, density = get_highest_density_kmer(G, seed_kmers[k], TMP_NODE_DENSITY)        
        if kmer is not None:
            best_kmer = kmer
            best_density = density
            break
    return best_kmer, best_density

def find_path(G, kmer_dict, kmax, seed):
    path = (seed,)
    density = G.node[seed][TMP_NODE_DENSITY]
    print "SEED", seed, density
    print "FORWARD"
    fwd_kmer, fwd_density = extend_seed(G, kmer_dict[FWD][path[-1]], kmax)
    while fwd_kmer is not None:
        path = path + fwd_kmer
        density = imin2(fwd_density, density)
        fwd_kmer, fwd_density = extend_seed(G, kmer_dict[FWD][path[-1]], kmax)
        print "NEW PATH", path
    print "REVERSE"
    print "KMERS", kmer_dict[REV][path[0]]
    rev_kmer, rev_density = extend_seed(G, kmer_dict[REV][path[0]], kmax)
    while rev_kmer is not None:
        path = rev_kmer + path
        density = imin2(rev_density, density)
        rev_kmer, rev_density = extend_seed(G, kmer_dict[REV][path[0]], kmax)
        print "NEW PATH", path
    print "FINAL PATH", path
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
        d[TMP_NODE_DENSITY] = imax2(MIN_DENSITY, d[TMP_NODE_DENSITY] - density)

def find_suboptimal_paths(G, kmer_dicts, 
                          kmax=3,
                          fraction_major_path=1e-3,
                          max_paths=100):
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
        path, density = find_path(G, kmer_dicts, kmax, seed)
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

