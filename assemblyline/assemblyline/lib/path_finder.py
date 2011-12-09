'''
Created on Dec 4, 2010

@author: mkiyer
'''
import logging
import collections
from assembler_base import NODE_DENSITY, EDGE_OUT_FRAC, EDGE_IN_FRAC, \
    MIN_NODE_DENSITY

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

def visit_node(G, src, dst):
    # get path attributes of src node
    src_attrs = G.node[src]
    src_length = src_attrs[PATH_LENGTH]
    src_weight = src_attrs[PATH_WEIGHT]
    # get dest node attributes
    dst_attrs = G.node[dst]
    dst_density = dst_attrs[TMP_NODE_DENSITY]
    dst_length = dst_attrs['length']
    dst_weight = dst_density * dst_length
    #assert dst_length == (dst.end - dst.start)
    # compute new path attributes when extending by this node
    length = src_length + dst_length
    # to compute the coverage going from the parent to the child, take
    # the existing path coverage at the parent node and multiply it by
    # the outgoing fraction, then add the child coverage multiplied by
    # its incoming fraction.
    edata = G.edge[src][dst]
    weight = (src_weight * edata[TMP_EDGE_OUT_FRAC] + 
              dst_weight * edata[TMP_EDGE_IN_FRAC])
    return weight, length

def dynprog_search(G, source):
    """
    Find the highest scoring path by dynamic programming    
    # Taken from NetworkX source code http://networkx.lanl.gov    
    """
    stack = collections.deque([source])
    G.node[source][PATH_NODES] = (source,)
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
            #logging.debug('\tweight=%f length=%d density=%f dst_density=%f' %
            #              (weight, length, density, dst_attrs[PATH_DENSITY]))
            # only update this node if it has not been visited
            # or if the new score is higher than the old score
            if ((dst_attrs[PATH_DENSITY] is None) or
                (density > dst_attrs[PATH_DENSITY])):
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

def subtract_path(G, path, density):
    """
    subtract the coverage density incurred by traversing the given path
    """
    imax2 = lambda x,y: x if x>=y else y
    for u in path:
        # subtract density but do not allow to become zero because that
        # will break some arithmetic steps
        d = G.node[u]
        d[TMP_NODE_DENSITY] = imax2(MIN_NODE_DENSITY, 
                                    d[TMP_NODE_DENSITY] - density) 

def find_suboptimal_paths(G, source, sink, fraction_major_path, 
                          max_paths=1000):
    """
    finds suboptimal paths through graph G using greedy algorithm
    that finds the highest density path using dynamic programming,
    removes that path, and then recomputes the best path. algorithm 
    stops after no more paths with positive density can be found,
    or when 'max_paths' paths have been generated
    """
    # copy the node weights to a temporary variable so that we can
    # manipulation them in the algorithm
    for n,d in G.nodes_iter(data=True):
        d[TMP_NODE_DENSITY] = d[NODE_DENSITY]
    for u,v,d in G.edges_iter(data=True):
        d[TMP_EDGE_OUT_FRAC] = d[EDGE_OUT_FRAC]
        d[TMP_EDGE_IN_FRAC] = d[EDGE_IN_FRAC]
    # find the best path
    logging.debug("\tFinding overall best path")
    path_results = {}
    init_path_attributes(G)
    dynprog_search(G, source)    
    # get resulting path and its attributes
    sink_attrs = G.node[sink]
    path = sink_attrs[PATH_NODES]
    weight = sink_attrs[PATH_WEIGHT]
    length = sink_attrs[PATH_LENGTH]
    density = sink_attrs[PATH_DENSITY]
    # unfortunately due to numerical issues the algorithm can 
    # rarely generate redundant paths.  to avoid this, we 
    # check each path against previous paths and ignore 
    # redundancies
    path_results[path] = density
    logging.debug("\t\tweight=%f length=%d density=%f" % 
                  (weight, length, density))
    # find suboptimal paths
    logging.debug("\tFinding suboptimal paths")
    density_lower_bound = fraction_major_path * density
    iterations = 0
    while (density >= density_lower_bound) and (iterations < max_paths):
        # remove influence of the "best" path
        subtract_path(G, path, density)
        # recompute edge fractions
        calculate_edge_fractions(G, TMP_NODE_DENSITY)
        # rerun dynamic programming to find next best path
        init_path_attributes(G)
        dynprog_search(G, source)
        # get path and store it
        sink_attrs = G.node[sink]
        path = sink_attrs[PATH_NODES]
        weight = sink_attrs[PATH_WEIGHT]
        length = sink_attrs[PATH_LENGTH]
        density = sink_attrs[PATH_DENSITY]
        logging.debug("\t\tweight=%f length=%d density=%f" % 
                      (weight, length, density))
        if path not in path_results:
            path_results[path] = density
        # count loops
        iterations += 1
    logging.debug("\t\tpath finding iterations=%d" % iterations)
    # cleanup graph attributes
    clear_path_attributes(G) 
    for n,d in G.nodes_iter(data=True):
        del d[TMP_NODE_DENSITY]
    for u,v,d in G.edges_iter(data=True):
        del d[TMP_EDGE_OUT_FRAC]
        del d[TMP_EDGE_IN_FRAC]
    # return (path,density) tuples
    return path_results.items()
