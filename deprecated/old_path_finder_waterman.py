'''
Created on Mar 15, 2011

@author: mkiyer
'''
'''
Created on Dec 4, 2010

@author: mkiyer
'''
import logging
import collections
import bisect
import networkx as nx

from base import NODE_WEIGHT, EDGE_OUT_FRAC, EDGE_IN_FRAC 

PLEN = 'plen'
PWEIGHT = 'pweight'
PDENSITY = 'pdensity'
PSRC = 'psrc'
PNPATHS = 'pnpaths'

def visit_node_reverse(G, src, dst):
    # compute length, coverage, and rpkm of the current
    # path when extended to include the child
    dst_attrs = G.node[dst]
    # new path length is parent length + child length
    length = dst_attrs.get(PLEN, 0) + (src.end - src.start)
    # to compute the coverage going from the parent to the child, take
    # the existing path coverage at the parent node and multiply it by
    # the outgoing fraction, then add the child coverage multiplied by
    # its incoming fraction.
    src_weight = G.node[src].get(NODE_WEIGHT, 0)
    dst_weight = dst_attrs.get(PWEIGHT, dst_attrs[NODE_WEIGHT])
    edata = G.edge[src][dst]
    weight = (src_weight * edata[EDGE_OUT_FRAC] + 
              dst_weight * edata[EDGE_IN_FRAC])
    return weight, length

def bfs_reverse(G, sink):
    """Compute the highest scoring path from each node to the
    end node

    Adapted from NetworkX source code http://networkx.lanl.gov
    """
    stack = collections.OrderedDict([(sink,1)])
    while stack:
        parent,parent_numpaths = stack.popitem()
        G.node[parent][PNPATHS] = parent_numpaths
        # iterate through children to find the maximum density
        # path from each child node to the sink
        for pred in G.predecessors_iter(parent):
            #logging.debug('Visiting %s <- %s' % (pred, parent))      
            weight, length = visit_node_reverse(G, pred, parent)
            # density is total coverage divided by the length of the path
            density = weight / float(length)
            #logging.debug('\twt=%f len=%d density=%f pred_density=%f' %
            #              (weight, length, density, G.node[pred].get(PDENSITY,-1)))
            # only update this node if it has not been visited
            # or if the new score is higher than the old score
            pred_attrs = G.node[pred]
            if density > pred_attrs.get(PDENSITY, -1):
                # keep pointer to parent node that produced this high scoring path
                pred_attrs[PLEN] = length
                pred_attrs[PWEIGHT] = weight
                pred_attrs[PDENSITY] = density
                pred_attrs[PSRC] = parent
                #logging.debug('\tupdated pred len=%d cov=%f density=%f' % 
                #              (length, weight, density))
            # continue iterating through the predecessor nodes
            if pred in stack:
                numpaths = stack[pred]
                del stack[pred]
            else:
                numpaths = parent_numpaths
            stack[pred] = 


def bfs_reverse(G, sink):
    """Compute the highest scoring path from each node to the
    end node

    Adapted from NetworkX source code http://networkx.lanl.gov
    """
    stack = [sink]
    # compute the number of paths from this node to the sink
    while stack:
        parent = stack[0]
        # compute number of paths from parent to sink
        npaths = sum(G.node[succ][PNPATHS] for succ in G.successors_iter(parent))
        if npaths == 0: npaths = 1
        G.node[parent][PNPATHS] = 1
        # iterate through children to find the maximum density
        # path from each child node to the sink
        for pred in G.predecessors_iter(parent):
            logging.debug('Visiting %s <- %s' % (pred, parent))      
            weight, length = visit_node_reverse(G, pred, parent)
            # density is total coverage divided by the length of the path
            density = weight / float(length)
            logging.debug('\twt=%f len=%d density=%f pred_density=%f' %
                          (weight, length, density, G.node[pred].get(PDENSITY,-1)))
            # only update this node if it has not been visited
            # or if the new score is higher than the old score
            pred_attrs = G.node[pred]
            if density > pred_attrs.get(PDENSITY, -1):
                # keep pointer to parent node that produced this high scoring path
                pred_attrs[PLEN] = length
                pred_attrs[PWEIGHT] = weight
                pred_attrs[PDENSITY] = density
                pred_attrs[PSRC] = parent
                logging.debug('\tupdated pred len=%d cov=%f density=%f' % 
                              (length, weight, density))
            # continue iterating through the predecessor node
            stack.append(pred)
        stack.pop(0)

def visit_node_forward(G, src, dst):
    src_attrs = G.node[src]
    length = src_attrs.get(PLEN, 0) + (dst.end - dst.start)
    src_weight = src_attrs.get(PWEIGHT, src_attrs[NODE_WEIGHT])
    dst_weight = G.node[dst].get(NODE_WEIGHT, 0)
    edata = G.edge[src][dst]
    weight = (src_weight * edata[EDGE_OUT_FRAC] + 
              dst_weight * edata[EDGE_IN_FRAC])
    return weight, length

def clear_path_attributes(G):
    '''
    remove dynamic programming node attributes that are 
    added to the graph temporarily
    
    must call this before calling the dynamic programming
    algorithm again on the same graph    
    '''
    # clear path attributes
    for n in G.nodes_iter():
        nattrs = G.node[n]
        if PLEN in nattrs:
            del nattrs[PLEN]
            del nattrs[PWEIGHT]
            del nattrs[PDENSITY]
            del nattrs[PSRC]
            del nattrs[PNPATHS]

def find_suboptimal_paths(G, source, sink, fraction):
    '''
    use dynamic programming to find high scoring suboptimal paths
    within some `fraction` of the best path
    '''
    # run the breath-first search to find the highest density 
    # from each node to the end node
    bfs_reverse(G, sink)
    # the score of the best scoring path is now stored at the source
    # node
    pattrs = G.node[source]
    density = pattrs[PDENSITY]
    length = pattrs[PLEN]
    weight = pattrs[PWEIGHT]
    # we will find suboptimal paths with density > fraction*best_density
    density_upper_bound = density
    density_lower_bound = density * fraction
    logging.debug("DENSITY upper=%f lower=%f" % 
                  (density_upper_bound, density_lower_bound))
    # depth first search to find all suboptimal paths
    stack = [([], source, 0, 0)]
    while stack:
        path, parent, pweight, plength = stack.pop()
        path.append(parent)
        if parent == sink:
            # reached the end of the graph
            logging.debug("RETURN path=%s" % (path))
            yield path, pweight, plength
            continue
        for child in G.successors_iter(parent):
            logging.debug('Visiting %s -> %s' % (parent, child))      
            # get child node length and weight
            cattrs = G.node[child]
            clength = (child.end - child.start)            
            cweight = cattrs.get(NODE_WEIGHT, 0)
            # compute density moving along this edge
            edata = G.edge[parent][child]
            t_weight = (pweight * edata[EDGE_OUT_FRAC] + 
                        cweight * edata[EDGE_IN_FRAC])
            t_length = plength + clength
            # get the upper bound of density achievable from child node
            # to the sink
            f_weight = cattrs.get(PWEIGHT, 0)
            f_length = cattrs.get(PLEN, 0)
            # now compute the prospective path density
            density = (t_weight + f_weight) / float(t_length + f_length)
            logging.debug('\tt(x,y)=%f,%f(%f) f(child)=%f,%f(%f)' % 
                          (t_weight, t_length, t_weight/t_length,
                           f_weight, f_length, f_weight/f_length))
            # if this density is greater than the lower bound, then the
            # current path remains viable for producing a suboptimal 
            # path
            if density >= density_lower_bound:
                # push this to stack
                stack.append((list(path), child, t_weight, t_length))
                logging.debug('\tadded child')
    # remove graph attributes
    clear_path_attributes(G)                    
