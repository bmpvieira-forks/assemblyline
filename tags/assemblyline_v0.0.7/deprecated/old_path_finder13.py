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
from base import imax2

# for dynamic programming algorithm
TMP_NODE_DENSITY = 'tmpnw'
TMP_EDGE_OUT_FRAC = 'tmpoutfrac'
TMP_EDGE_IN_FRAC = 'tmpinfrac'
PATH_DICT = 'ppdict'

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
    for u,v,d in G.edges_iter(data=True):
        d[TMP_EDGE_OUT_FRAC] = d[EDGE_OUT_FRAC]
        d[TMP_EDGE_IN_FRAC] = d[EDGE_IN_FRAC]

def clear_tmp_attrs(G):
    for n,d in G.nodes_iter(data=True):
        del d[TMP_NODE_DENSITY]
    for u,v,d in G.edges_iter(data=True):
        del d[TMP_EDGE_OUT_FRAC]
        del d[TMP_EDGE_IN_FRAC]

def init_path_attributes(G):
    """
    must call this before calling the dynamic programming
    algorithm
    """
    # reset path attributes
    for n,d in G.nodes_iter(data=True):
        d[PATH_DICT] = {}

def clear_path_attributes(G):
    '''
    remove dynamic programming node attributes that are 
    added to the graph temporarily
    
    must call this before calling the dynamic programming
    algorithm again on the same graph    
    '''
    # clear path attributes
    for n,d in G.nodes_iter(data=True):
        del d[PATH_DICT]

def visit_node(G, v):
    # get current node attributes
    v_attrs = G.node[v]
    v_length = v_attrs[NODE_LENGTH]
    v_weight = v_length * v_attrs[TMP_NODE_DENSITY]
    for w in G.successors_iter(v):
        vw_out_frac = G[v][w][TMP_EDGE_OUT_FRAC]
        if len(v_attrs[PATH_DICT]) == 0:
            # node has no predecessors
            best_path_weight = v_weight * vw_out_frac
            best_path_length = v_length
        else:
            # iterate through predecessors to find the best predecessor
            best_path_density = None
            for u,path_info in v_attrs[PATH_DICT].iteritems():
                uv_in_frac = G[u][v][TMP_EDGE_IN_FRAC]
                path_weight, path_length = path_info
                # calculate weight and density of node 'v' given u -> v -> w path
                new_path_weight = path_weight + (uv_in_frac * v_weight * vw_out_frac)
                new_path_length = path_length + v_length
                new_path_density = new_path_weight / float(new_path_length) 
                # update maximum path density seen when passing u -> v -> w
                if (best_path_density is None) or (new_path_density > best_path_density):
                    best_path_weight = new_path_weight
                    best_path_length = float(new_path_length)
                    best_path_density = best_path_weight / float(best_path_length)
        # store best path density for v -> w in 'w' dictionary
        w_attrs = G.node[w]
        w_attrs[PATH_DICT][v] = (best_path_weight, best_path_length)

def dynprog_search(G, source, reverse=False):
    """
    Find the highest scoring path by dynamic programming    
    # Taken from NetworkX source code http://networkx.lanl.gov    
    """
    # setup initial path attributes
    init_path_attributes(G)
    queue = collections.deque(nx.topological_sort(G))
    # only consider nodes reachable from source
    path_length_dict = nx.shortest_path_length(G, source=source)
    while queue:
        current = queue.popleft()
        # ensure there is a path from source to node
        if current not in path_length_dict:
            continue
        visit_node(G, current)

def get_best_pred(pred_dict):
    max_pred = None
    max_density = None
    max_weight = None
    max_length = None
    for pred,wtlen in pred_dict.iteritems():
        weight, length = wtlen
        if length == 0:
            density = 0.0
        else:
            density = weight / float(length)
        if (max_density is None) or (density > max_density):
            max_pred = pred
            max_density = density
            max_weight = weight
            max_length = length
    return max_pred, max_weight, max_length

def traceback(G, sink, inclusive=True):
    """
    compute path weight and length, and subtract density
    from nodes along path 
    """
    # get predecessor node with max density
    pred_dict = G.node[sink][PATH_DICT]
    pred, path_weight, path_length = get_best_pred(pred_dict)
    path = []
    if inclusive:
        # add sink node weight, length
        path.append(sink)
        sink_attrs = G.node[sink]
        sink_length = sink_attrs[NODE_LENGTH]
        sink_weight = sink_length * sink_attrs[TMP_NODE_DENSITY]
        sink_weight *= G[pred][sink][TMP_EDGE_IN_FRAC]
        path_weight += sink_weight
        path_length += sink_length
    # iterate over predecessor nodes
    while len(pred_dict) > 0:
        # get predecessor node with max density
        pred, weight, length = get_best_pred(pred_dict)
        path.append(pred)
        pred_dict = G.node[path[-1]][PATH_DICT]
    # path is now discovered in reverse order
    path.reverse()
    return tuple(path), path_weight, path_length

def reverse_graph(G):
    # reverse all graph edges so that search goes in
    # opposite direction
    G.reverse(copy=False)
    # flip in/out edge fractions
    for u,nbrsdict in G.adjacency_iter():
        for v,eattr in nbrsdict.items():
            tmp = eattr[TMP_EDGE_OUT_FRAC]                
            eattr[TMP_EDGE_OUT_FRAC] = eattr[TMP_EDGE_IN_FRAC] 
            eattr[TMP_EDGE_IN_FRAC] = tmp

def find_path(G, source, sink, reverse=False, inclusive=True):
    """
    G - graph
    source, sink - start/end nodes
    """    
    if reverse:
        source, sink = sink, source
        reverse_graph(G)
    # dynamic programming search for best path
    dynprog_search(G, source)
    # traceback to get path
    path, weight, length = traceback(G, sink, inclusive)
    # restore edges
    if reverse:
        source, sink = sink, source
        reverse_graph(G)
    return path, weight, length

def get_seed_subgraph(H, A, nodes, node_indexes, seed):
    seed_index = node_indexes[seed]
    subgraph_indexes = (A[:,seed_index].getA1() + 
                        A[seed_index,:].getA1()).nonzero()[0]
    subgraph_nodes = tuple(nodes[i] for i in subgraph_indexes)
    G = H.subgraph(subgraph_nodes)
    return G

def subtract_path(G, path):
    """
    subtract density from nodes along path 
    """
    # first and last node are 'dummy' nodes so don't need to
    # include them in the loop
    for i in xrange(1, len(path)-1):        
        u = path[i-1]
        v = path[i]
        w = path[i+1]        
        d = G.node[v]
        # compute length and weight up to this point along path
        v_density_used = (G[u][v][TMP_EDGE_IN_FRAC] * 
                          G[v][w][TMP_EDGE_OUT_FRAC] * 
                          d[TMP_NODE_DENSITY])
        # subtract density
        d[TMP_NODE_DENSITY] = imax2(MIN_NODE_DENSITY,
                                    d[TMP_NODE_DENSITY] - v_density_used)

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
    # find suboptimal paths
    logging.debug("\tFinding best paths in graph with %d nodes" % (len(G)))
    path, weight, length = find_path(G, source, sink)
    density = weight / float(length)
    lowest_density = density * fraction_major_path
    logging.debug("\t\tweight=%f length=%d density=%f" % 
                  (weight, length, density))
    # subtract density of path
    subtract_path(G, path)
    # recompute edge fractions
    calculate_edge_fractions(G, TMP_NODE_DENSITY, 
                             in_frac_attr=TMP_EDGE_IN_FRAC,
                             out_frac_attr=TMP_EDGE_OUT_FRAC)
    # store path
    if path not in path_results:
        path_results[path] = density
    logging.debug("\tFinding suboptimal paths in graph with %d nodes" % 
                  (len(G)))
    iterations = 0    
    while iterations < max_paths:
        # find next path and update graph
        path, weight, length = find_path(G, source, sink)
        density = weight / float(length)
        # break out of loop if density is too low
        if density < lowest_density:
            break
        logging.debug("\t\tweight=%f length=%d density=%f" % 
                      (weight, length, density))
        # subtract density of path
        subtract_path(G, path)
        # recompute edge fractions
        calculate_edge_fractions(G, TMP_NODE_DENSITY, 
                                 in_frac_attr=TMP_EDGE_IN_FRAC,
                                 out_frac_attr=TMP_EDGE_OUT_FRAC)
        # store path
        if path not in path_results:
            path_results[path] = density
        iterations +=1
    logging.debug("\t\tpath finding iterations=%d" % iterations)
    # cleanup graph attributes
    clear_path_attributes(G) 
    clear_tmp_attrs(G)
    # return (path,density) tuples sorted from high -> low density
    return sorted(path_results.items(), key=operator.itemgetter(1), reverse=True)
    #path_results = [(p,d) for (p,d) in path_results.items() if 
    #                d >= lowest_density]
    #path_results.sort(key=operator.itemgetter(1), reverse=True)
    #return path_results
