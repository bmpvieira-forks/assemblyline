'''
Created on Dec 4, 2010

@author: mkiyer
'''
import logging
import collections
import networkx as nx

from assembler_base import NODE_WEIGHT, EDGE_OUT_FRAC, EDGE_IN_FRAC, EDGE_DENSITY

# for get num paths algorithm
PNPATHS = 'pnpaths'

# for dynamic programming algorithm
TMP_NODE_WEIGHT = 'tmpnw'
TMP_EDGE_DENSITY = 'tmpdensity'
TMP_EDGE_OUT_FRAC = 'tmpoutfrac'
TMP_EDGE_IN_FRAC = 'tmpinfrac'
PATH_LENGTH = 'plen'
PATH_WEIGHT = 'pweight'
PATH_DENSITY = 'pdensity'
PATH_NODES = 'pnodes'

# constant minimum edge density
MIN_EDGE_DENSITY = 0.001

def visit_node(G, src, dst):
    # get path attributes of src node
    src_attrs = G.node[src]
    src_length = src_attrs[PATH_LENGTH]
    src_weight = src_attrs[PATH_WEIGHT]
    # get dest node attributes
    dst_weight = G.node[dst][TMP_NODE_WEIGHT]
    dst_length = (dst.end - dst.start)
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

def traceback(G, source, sink):
    """
    perform traceback to find the highest density path
    """
    sink_attrs = G.node[sink]
    path = sink_attrs[PATH_NODES]
    tmp_length = sink_attrs[PATH_LENGTH]
    tmp_weight = sink_attrs[PATH_WEIGHT]
    # reconstitute path, length, and weight
    weight = 0.0
    length = 0
    for i in xrange(1, len(path)):
        src = path[i-1]
        dst = path[i]
        dst_weight = G.node[dst][NODE_WEIGHT]
        # compute length and weight
        length += (dst.end - dst.start)
        edata = G.edge[src][dst]
        weight = (weight * edata[EDGE_OUT_FRAC] + 
                  dst_weight * edata[EDGE_IN_FRAC])
    #logging.debug("TRACEBACK path=%s numnodes=%d tmp_weight=%f tmp_length=%d weight=%f length=%d" % 
    #              (path, len(path), tmp_weight, tmp_length, weight, length))
    assert tmp_length == length
    return path, tmp_weight, weight, length

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

def recalculate_edge_attrs(G, u, v):
    # find total weight leaving node 'u'
    u_succs = G.successors(u)
    out_density_total = sum(G.edge[u][x][TMP_EDGE_DENSITY] for x in u_succs)
    assert out_density_total > 0
    # compute edge attributes leaving node 'u'
    for succ in u_succs:
        eattrs = G.edge[u][succ]
        out_frac = eattrs[TMP_EDGE_DENSITY]/out_density_total
        eattrs[TMP_EDGE_OUT_FRAC] = out_frac
    # find total density entering node 'v'
    v_preds = G.predecessors(v)
    in_density_total = sum(G.edge[x][v][TMP_EDGE_DENSITY] for x in v_preds)
    assert in_density_total > 0
    # compute edge density fractions entering node 'v'
    for pred in v_preds:
        eattrs = G.edge[pred][v]
        in_frac = eattrs[TMP_EDGE_DENSITY]/in_density_total
        eattrs[TMP_EDGE_IN_FRAC] = in_frac

def subtract_path(G, path, weight, length):
    """
    subtract the coverage density incurred by traversing the given path
    """
    # find length of nonzero nodes in path
    pos_length = 0
    pos_weight = 0
    for v in path:
        if G.node[v][TMP_NODE_WEIGHT] > 0:
            pos_length += (v.end - v.start)
            pos_weight += G.node[v][TMP_NODE_WEIGHT] 
    if pos_length == 0:
        return 0
    density = max(weight, pos_weight) / float(pos_length)   
    for i,v in enumerate(path):
        v_weight_used = density * (v.end - v.start)
        new_weight = G.node[v][TMP_NODE_WEIGHT] - v_weight_used        
        G.node[v][TMP_NODE_WEIGHT] = max(0, new_weight)
        #G.node[v][TMP_NODE_WEIGHT] -= v_weight_used
        if i == 0:
            continue
        # subtract density from edge
        u = path[i-1]
        eattrs = G.edge[u][v]
        new_edge_density = eattrs[TMP_EDGE_DENSITY] - density
        eattrs[TMP_EDGE_DENSITY] = max(MIN_EDGE_DENSITY, new_edge_density)
        # recompute edge attrs
        recalculate_edge_attrs(G, u, v)
    return pos_weight

def get_num_paths(G, source, sink):
    G.reverse(copy=False)
    ordered_nodes = nx.topological_sort(G)
    G.node[ordered_nodes[0]][PNPATHS] = 1
    for n in ordered_nodes[1:]:
        G.node[n][PNPATHS] = sum(G.node[pred][PNPATHS] for pred in G.predecessors_iter(n))
    G.reverse(copy=False)

def find_suboptimal_paths(G, source, sink, 
                          max_exhaustive=1000):
    """
    use dynamic programming to find high scoring suboptimal paths
    within some `fraction` of the best path
    """
    # compute the number of total paths in the graph
    #get_num_paths(G, source, sink)
    #pattrs = G.node[source]
    #numpaths = pattrs[PNPATHS]
    #logging.debug("\ttotal_paths=%d" % (numpaths))
    # copy the node weights to a temporary variable so that we can
    # manipulation them in the algorithm
    for n,d in G.nodes_iter(data=True):
        d[TMP_NODE_WEIGHT] = d[NODE_WEIGHT]
    for u,v,d in G.edges_iter(data=True):
        d[TMP_EDGE_DENSITY] = d[EDGE_DENSITY]
        d[TMP_EDGE_OUT_FRAC] = d[EDGE_OUT_FRAC]
        d[TMP_EDGE_IN_FRAC] = d[EDGE_IN_FRAC]
    # find the best path
    logging.debug("\tFinding overall best path")
    path_results = {}
    init_path_attributes(G)
    dynprog_search(G, source)    
    path, tmp_weight, weight, length = traceback(G, source, sink)
    path_results[path] = (weight, length)
    logging.debug("\t\ttmp_weight=%f weight=%f length=%d density=%f" % 
                  (tmp_weight, weight, length, weight/float(length)))
    # find suboptimal paths
    logging.debug("\tFinding suboptimal paths")
    iterations = 0
    while tmp_weight > 0:
        # remove influence of the "best" path
        subtract_path(G, path, weight, length)
        # rerun dynamic programming to find next best path
        init_path_attributes(G)
        dynprog_search(G, source)
        path, tmp_weight, weight, length = traceback(G, source, sink)
        if path not in path_results:
            path_results[path] = (weight, length)
        logging.debug("\t\ttmp_weight=%f weight=%f length=%d density=%f" % 
                      (tmp_weight, weight, length, weight/float(length)))
        # count loops
        iterations += 1
    logging.debug("\t\tpath finding iterations=%d" % iterations)
    for path_tuple, res in path_results.iteritems():
        weight, length = res
        yield path_tuple, weight, length
    # cleanup graph attributes
    clear_path_attributes(G) 
    for n,d in G.nodes_iter(data=True):
        del d[TMP_NODE_WEIGHT]
    for u,v,d in G.edges_iter(data=True):
        del d[TMP_EDGE_DENSITY]
        del d[TMP_EDGE_OUT_FRAC]
        del d[TMP_EDGE_IN_FRAC]
