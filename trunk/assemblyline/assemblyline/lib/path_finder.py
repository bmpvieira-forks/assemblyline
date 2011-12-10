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
PATH_LENGTH = 'plen'
PATH_WEIGHT = 'pweight'
PATH_DENSITY = 'pdensity'
PATH_NODES = 'pnodes'

def init_path_attributes(G):
    """
    must call this before calling the dynamic programming
    algorithm
    """
    # reset path attributes
    for n,d in G.nodes_iter(data=True):
        d[PATH_LENGTH] = 0
        d[PATH_WEIGHT] = 0
        d[PATH_DENSITY] = None
        d[PATH_NODES] = tuple()

def clear_path_attributes(G):
    '''
    remove dynamic programming node attributes that are 
    added to the graph temporarily
    
    must call this before calling the dynamic programming
    algorithm again on the same graph    
    '''
    # clear path attributes
    for n,d in G.nodes_iter(data=True):
        del d[PATH_LENGTH]
        del d[PATH_WEIGHT]
        del d[PATH_DENSITY]
        del d[PATH_NODES]

def visit_node(G, src, dst, reverse=False):
    # get path attributes of src node
    src_attrs = G.node[src]
    src_length = src_attrs[PATH_LENGTH]
    src_weight = src_attrs[PATH_WEIGHT]
    # get dest node attributes
    dst_attrs = G.node[dst]
    dst_density = dst_attrs[TMP_NODE_DENSITY]
    dst_length = dst_attrs[NODE_LENGTH]
    dst_weight = dst_density * dst_length
    # compute new path attributes when extending by this node
    length = src_length + dst_length
    # to compute the coverage going from the parent to the child, take
    # the existing path coverage at the parent node and multiply it by
    # the outgoing fraction, then add the child coverage multiplied by
    # its incoming fraction.
    edata = G[src][dst]
    if reverse:
        in_frac, out_frac = edata[TMP_EDGE_OUT_FRAC], edata[TMP_EDGE_IN_FRAC]
    else:
        in_frac, out_frac = edata[TMP_EDGE_IN_FRAC], edata[TMP_EDGE_OUT_FRAC]
    weight = (src_weight * out_frac +
              dst_weight * in_frac)
    return weight, length

def dynprog_search(G, source, reverse=False):
    """
    Find the highest scoring path by dynamic programming    
    # Taken from NetworkX source code http://networkx.lanl.gov    
    """
    # initialize node attributes
    init_path_attributes(G)
    # add source attributes
    stack = collections.deque([source])
    src_attrs = G.node[source]
    src_attrs[PATH_LENGTH] = src_attrs[NODE_LENGTH]
    src_attrs[PATH_WEIGHT] = src_attrs[NODE_LENGTH] * src_attrs[NODE_DENSITY]
    src_attrs[PATH_DENSITY] = src_attrs[NODE_DENSITY]
    src_attrs[PATH_NODES] = (source,)
    if reverse:
        # reverse all graph edges so that search goes in
        # opposite direction
        G.reverse(copy=False)    
    while stack:
        current = stack[0]
        current_attrs = G.node[current]
        for dst in G.successors_iter(current):
            # compute weight, length, and density at dest when visited
            # from the current node
            #logging.debug('Visiting %s -> %s' % (current, dst))            
            weight, length = visit_node(G, current, dst)
            density = weight / float(length)
            dst_attrs = G.node[dst]
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
        G.reverse(copy=False)

#def traceback(G, source, sink):
#    """
#    perform traceback to find the highest density path
#    """
#    sink_attrs = G.node[sink]
#    path = sink_attrs[PATH_NODES]
#    # reconstitute path, length, and weight
#    weight = 0.0
#    length = 0
#    for i in xrange(1, len(path)):
#        src = path[i-1]
#        # get dst node attributes
#        dst = path[i]
#        dst_attrs = G.node[dst]
#        dst_density = dst_attrs[NODE_DENSITY]
#        dst_length = dst_attrs['length']
#        dst_weight = dst_density * float(dst_length)
#        # compute length and weight up to this point along path
#        length += dst_length
#        edata = G.edge[src][dst]
#        weight = (weight * edata[EDGE_OUT_FRAC] + 
#                  dst_weight * edata[EDGE_IN_FRAC])
#    return path, weight, length

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

#def subtract_path(G, path, density):
#    """
#    subtract the coverage density incurred by traversing the given path
#    """
#    for u in path:
#        # subtract density but do not allow to become zero because that
#        # will break some arithmetic steps
#        d = G.node[u]
#        d[TMP_NODE_DENSITY] = imax2(MIN_NODE_DENSITY, 
#                                    d[TMP_NODE_DENSITY] - density) 

def subtract_path(G, path):
    """
    subtract density from each node implied by the path
    """
    for i in xrange(1, len(path)-1):        
        u = path[i-1]
        v = path[i]
        w = path[i+1]        
        d = G.node[v]
        density_used = (G[u][v][EDGE_IN_FRAC] * 
                        G[v][w][EDGE_OUT_FRAC] * 
                        d[TMP_NODE_DENSITY])
        d[TMP_NODE_DENSITY] = imax2(MIN_NODE_DENSITY,
                                    d[TMP_NODE_DENSITY] - density_used)

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
    fwd_weight = sink_attrs[PATH_WEIGHT]
    fwd_length = sink_attrs[PATH_LENGTH]
    # get best path from source to node
    #logging.debug("\tsearching reverse")
    dynprog_search(G, start_node, reverse=True)
    source_attrs = G.node[source]
    rev_path = source_attrs[PATH_NODES]
    rev_weight = source_attrs[PATH_WEIGHT]
    rev_length = source_attrs[PATH_LENGTH]
    # merge paths
    assert rev_path[0] == fwd_path[0]
    path = rev_path[::-1] + fwd_path[1:]
    start_node_data = G.node[start_node]
    start_node_length = start_node_data[NODE_LENGTH]
    start_node_weight = start_node_length * start_node_data[TMP_NODE_DENSITY]
    length = fwd_length + rev_length - start_node_length
    weight = fwd_weight + rev_weight - start_node_weight
    return path, length, weight

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
    node_data_list = G.nodes()
    node_data_list.remove(source)
    node_data_list.remove(sink)
    node_data_list = [(n,G.node[n]) for n in node_data_list]
    logging.debug("\tFinding suboptimal paths in graph with %d nodes" % 
                  (len(node_data_list)))
    iterations = 0
    highest_density = 0.0
    while (len(node_data_list) > 0) and (iterations < max_paths):
        # maintain list of node-data tuples in sorted order
        # such that last item in list is highest density node
        node_data_list.sort(key=lambda nd: nd[1][TMP_NODE_DENSITY])
        # pop the "best" node from list and use it as a seed for path
        best_node, best_node_data = node_data_list.pop()
        if best_node_data[TMP_NODE_DENSITY] < (highest_density * fraction_major_path):
            break
        logging.debug("\t\thighest density node %s density=%f" % 
                      (best_node, best_node_data[TMP_NODE_DENSITY]))
        # dynamic programming search for best path that includes this node
        path, length, weight = find_path(G, best_node, source, sink)
        density = weight / float(length)
        logging.debug("\t\tweight=%f length=%d density=%f" % 
                      (weight, length, density))
        # store path
        if path not in path_results:
            path_results[path] = density
        # update highest density
        if density > highest_density:
            highest_density = density
        # remove the density of the path from the graph
        subtract_path(G, path)
        # recompute edge fractions
        calculate_edge_fractions(G, TMP_NODE_DENSITY)
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