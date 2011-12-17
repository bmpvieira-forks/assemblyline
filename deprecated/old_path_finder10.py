'''
Created on Dec 9, 2011

@author: mkiyer
'''
import logging
import collections
import operator
from assembler_base import NODE_DENSITY, NODE_LENGTH, EDGE_OUT_FRAC, EDGE_IN_FRAC, \
    MIN_NODE_DENSITY
from base import imax2

# for dynamic programming algorithm
TMP_NODE_DENSITY = 'tmpnw'
TMP_EDGE_OUT_FRAC = 'tmpoutfrac'
TMP_EDGE_IN_FRAC = 'tmpinfrac'
PATH_DICT = 'ppdict'

def calculate_edge_fractions(G, node_density_attr=NODE_DENSITY):
    for u,nbrdict in G.adjacency_iter():
        # find total coverage flowing out of this node
        out_density_total = sum(G.node[v][node_density_attr] for v in nbrdict)
        for v, eattrs in nbrdict.iteritems():
            # find fraction flowing out of each node
            out_frac = G.node[v][node_density_attr]/float(out_density_total)
            eattrs[EDGE_OUT_FRAC] = out_frac            
    # now fraction of weight flowing in to each node
    # reverse edge directions and use adjacency iter again
    G.reverse(copy=False)
    for u,nbrdict in G.adjacency_iter():
        # find total coverage density flowing into this node
        in_density_total = sum(G.node[v][node_density_attr] for v in nbrdict)
        for v, eattrs in nbrdict.iteritems():
            # find fraction flowing into node            
            in_frac = G.node[v][node_density_attr]/float(in_density_total)
            eattrs[EDGE_IN_FRAC] = in_frac
    # reverse the edges back to normal
    G.reverse(copy=False)

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

def dynprog_search(G, source, reverse=False):
    """
    Find the highest scoring path by dynamic programming    
    # Taken from NetworkX source code http://networkx.lanl.gov    
    """
    # initialize node attributes
    init_path_attributes(G)
    if reverse:
        reverse_graph(G)
    # setup initial path attributes
    stack = collections.deque([source])
    while stack:
        current = stack[0]
        visit_node(G, current)
        stack.extend(G.successors(current))
        stack.popleft()
        
        
        
        current_attrs = G.node[current]
        for dst in G.successors_iter(current):
            # compute weight, length, and density at dest when visited
            # from the current node
            #logging.debug('Visiting %s -> %s' % (current, dst))            
            weight, length = visit_node(G, current, dst)
            density = weight / float(length)
            # compare to dest node density
            dst_attrs = G.node[dst]
            if len(dst_attrs[PATH_NODES]) == 0:
                dst_density = None
            else:
                dst_density = dst_attrs[PATH_WEIGHT] / dst_attrs[PATH_LENGTH]
            #logging.debug('\tweight=%f length=%d density=%f' %
            #              (weight, length, density))
            # update this node if it has not been visited
            # or if the density is greater than previously seen.
            # break density ties by taking path with greater total
            # coverage  
            if ((dst_attrs[PATH_DENSITY] is None) or 
                (density > dst_attrs[PATH_DENSITY]) or
                ((density == dst_attrs[PATH_DENSITY]) and
                 (weight > dst_attrs[PATH_WEIGHT]))):
                # keep pointer to current node that produced this 
                # high scoring path
                dst_attrs[PATH_LENGTH] = length
                dst_attrs[PATH_WEIGHT] = weight
                dst_attrs[PATH_DENSITY] = density
                dst_attrs[PATH_NODES] = current_attrs[PATH_NODES] + (dst,)
                #logging.debug('\tupdated dst weight=%f length=%d density=%f' % 
                #              (weight, length, density))
                # continue iterating through the dst path
                stack.append(dst)
        stack.popleft()
    # restore edges
    if reverse:
        reverse_graph(G)

def traceback(G, path):
    """
    compute path weight and length, and subtract density
    from nodes along path 
    """
    # reconstitute path, length, and weight
    weight = 0.0
    length = 0
    # first and last node are 'dummy' nodes so don't need to
    # include them in the loop
    for i in xrange(1, len(path)-1):        
        u = path[i-1]
        v = path[i]
        w = path[i+1]        
        d = G.node[v]
        # compute length and weight up to this point along path
        v_length = d[NODE_LENGTH]        
        v_density_used = (G[u][v][EDGE_IN_FRAC] * 
                          G[v][w][EDGE_OUT_FRAC] * 
                          d[TMP_NODE_DENSITY])
        # update total length and weight
        length += v_length
        weight += v_length * v_density_used
        # subtract density
        d[TMP_NODE_DENSITY] = imax2(MIN_NODE_DENSITY,
                                    d[TMP_NODE_DENSITY] - v_density_used)
    return weight, length

def find_path(G, start_node, source, sink):
    """
    find best path from source -> node -> sink
    """
    # get best path from node to sink
    #logging.debug("\tsearching forward")
    dynprog_search(G, start_node)
    # get resulting path and its attributes
    sink_attrs = G.node[sink]
    fwd_path = sink_attrs[PATH_NODES]
    # get best path from source to node
    #logging.debug("\tsearching reverse")
    dynprog_search(G, start_node, reverse=True)
    source_attrs = G.node[source]
    rev_path = source_attrs[PATH_NODES]
    # merge paths
    assert rev_path[0] == fwd_path[0]
    path = rev_path[::-1] + fwd_path[1:]
    return path


def find_suboptimal_paths(G, source, sink, fraction_major_path,
                          max_paths=100):
    # copy the node weights to a temporary variable so that we can
    # manipulation them in the algorithm
    for n,d in G.nodes_iter(data=True):
        d[TMP_NODE_DENSITY] = d[NODE_DENSITY]
    for u,v,d in G.edges_iter(data=True):
        d[TMP_EDGE_OUT_FRAC] = d[EDGE_OUT_FRAC]
        d[TMP_EDGE_IN_FRAC] = d[EDGE_IN_FRAC]


def find_suboptimal_paths(G, source, sink, fraction_major_path,
                          max_paths=100):
    """
    finds suboptimal paths through graph G using greedy algorithm
    that chooses the highest density node, extends a full path
    that includes that node and the source/sink, returns the path,
    subtracts the path density from the graph, and repeats
    
    algorithm stops when path density falls below a fraction of the
    highest density path found in the graph, or when 'max_paths'
    iterations have completed
    """
    # copy the node weights to a temporary variable so that we can
    # manipulation them in the algorithm
    for n,d in G.nodes_iter(data=True):
        d[TMP_NODE_DENSITY] = d[NODE_DENSITY]
    for u,v,d in G.edges_iter(data=True):
        d[TMP_EDGE_OUT_FRAC] = d[EDGE_OUT_FRAC]
        d[TMP_EDGE_IN_FRAC] = d[EDGE_IN_FRAC]
    # store paths in a dictionary in order to avoid redundant paths
    # that arise due to issues in the data
    path_results = collections.OrderedDict()
    # get starting list of nodes and remove source/sink
    used_nodes = set(G.nodes())
    used_nodes.remove(source)
    used_nodes.remove(sink)
    logging.debug("\tFinding suboptimal paths in graph with %d nodes" % 
                  (len(G)))
    iterations = 0
    highest_density = 0.0
    while (len(used_nodes) > 0) and (iterations < max_paths):
        # dynamic programming search for best path
        dynprog_search(G, source)
        # get path
        G.node[sink][PATH_NODES]

        # get path and store it
        sink_attrs = G.node[sink]
        path = sink_attrs[PATH_NODES]
        weight = sink_attrs[PATH_WEIGHT]
        length = sink_attrs[PATH_LENGTH]
        density = sink_attrs[PATH_DENSITY]
        if (density < density_lower_bound):
            break 
        logging.debug("\t\tweight=%f length=%d density=%f" % 
                      (weight, length, density))
        
        
        path = find_path(G, seed_node, source, sink)
        # process path to find weight and length, and update node weights
        weight, length = traceback(G, path)
        density = weight / float(length)
        logging.debug("\t\tweight=%f length=%d density=%f" % 
                      (weight, length, density))
        # recompute edge fractions
        calculate_edge_fractions(G, TMP_NODE_DENSITY)
        # store path
        if path not in path_results:
            path_results[path] = density
        # update highest density
        if density > highest_density:
            highest_density = density
        iterations +=1
    logging.debug("\t\tpath finding iterations=%d" % iterations)
    # cleanup graph attributes
    clear_path_attributes(G) 
    for n,d in G.nodes_iter(data=True):
        del d[TMP_NODE_DENSITY]
    for u,v,d in G.edges_iter(data=True):
        del d[TMP_EDGE_OUT_FRAC]
        del d[TMP_EDGE_IN_FRAC]
    # return (path,density) tuples sorted from high -> low density
    lowest_density = highest_density * fraction_major_path
    path_results = [(p,d) for (p,d) in path_results.items() if 
                    d >= lowest_density]
    path_results.sort(key=operator.itemgetter(1), reverse=True)
    return path_results