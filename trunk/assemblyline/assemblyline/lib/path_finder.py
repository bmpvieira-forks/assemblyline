'''
Created on Dec 4, 2010

@author: mkiyer
'''
import logging
import collections
import bisect
import random
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

def dynprog_reverse(G, sink):
    """Compute the highest scoring path from each node to the
    end node

    Adapted from NetworkX source code http://networkx.lanl.gov
    """
    stack = [sink]
    # compute the number of paths from this node to the sink
    while stack:
        parent = stack[0]
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
                # continue iterating through the predecessor node
                stack.append(pred)
        stack.pop(0)

def traceback(G, source, sink):
    """perform traceback to find the highest density path"""
    path = [source]
    pattrs = G.node[source]
    weight = pattrs[PWEIGHT]    
    length = pattrs[PLEN]
    while path[-1] != sink:
        #logging.debug('path=%s node=%s attrs=%s' %
        #              (path, path[-1], G.node[path[-1]]))
        path.append(G.node[path[-1]][PSRC])
    return path, weight, length

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

def get_num_paths(G, source, sink):
    G.reverse(copy=False)
    ordered_nodes = nx.topological_sort(G)
    G.node[ordered_nodes[0]][PNPATHS] = 1
    for n in ordered_nodes[1:]:
        G.node[n][PNPATHS] = sum(G.node[pred][PNPATHS] for pred in G.predecessors_iter(n))
    G.reverse(copy=False)

def find_paths_dfs(G, source, sink, max_paths):
    '''
    use depth-first search to exhaustively find all paths
    in the graph
    '''
    # depth first search to find all paths
    stack = [([], source, 0, 0)]
    numpaths = 0
    while stack:
        path, parent, pweight, plength = stack.pop()
        path.append(parent)
        if parent == sink:
            # reached the end of the graph
            #logging.debug("RETURN path=%s" % (path))
            yield path, pweight, plength
            numpaths += 1
            if numpaths > max_paths:
                break
            continue
        for child in G.successors_iter(parent):
            #logging.debug('Visiting %s -> %s' % (parent, child))      
            # get child node length and weight
            cattrs = G.node[child]
            clength = (child.end - child.start)            
            cweight = cattrs.get(NODE_WEIGHT, 0)
            # compute density moving along this edge
            edata = G.edge[parent][child]
            t_weight = (pweight * edata[EDGE_OUT_FRAC] + 
                        cweight * edata[EDGE_IN_FRAC])
            t_length = plength + clength
            # push this to stack
            stack.append((list(path), child, t_weight, t_length))
            #logging.debug('\tadded child')


def choose_successor(G, succs):
    succs_numpaths = [G.node[succ][PNPATHS] for succ in succs]
    total_numpaths = sum(succs_numpaths)
    r = random.randrange(total_numpaths)
    totalcount = 0
    for i,count in enumerate(succs_numpaths):
        totalcount += count
        if totalcount < r:
            break
    #logging.debug("Chose %d/%d successors" % (i+1, len(succs)))
    return succs[i]

def find_paths_markov(G, source, sink, maxpaths):
    """
    randomly sample N paths in the graph
    """
    # randomly sample paths with depth first search
    numpaths = 0
    while numpaths < maxpaths:
        path = [source]
        length = 0
        weight = 0
        while path[-1] != sink:
            # choose a random successor node to visit
            child = choose_successor(G, G.successors(path[-1]))
            # get child node length and weight
            cattrs = G.node[child]
            clength = (child.end - child.start)            
            cweight = cattrs.get(NODE_WEIGHT, 0)
            # compute weight and length moving along this edge
            edata = G.edge[path[-1]][child]
            weight = (weight * edata[EDGE_OUT_FRAC] + 
                      cweight * edata[EDGE_IN_FRAC])
            length += clength
            # add child
            path.append(child)
        yield path, weight, length
        numpaths += 1

def find_suboptimal_paths(G, source, sink, max_exhaustive=1000):
    '''
    use dynamic programming to find high scoring suboptimal paths
    within some `fraction` of the best path
    '''
    # compute the number of total paths in the graph
    logging.debug("Graph (%s -> %s)" % (source, sink))
    get_num_paths(G, source, sink)
    # decide whether to compute suboptimal paths using random
    # walking or by exhaustive search
    pattrs = G.node[source]
    numpaths = pattrs[PNPATHS]
    logging.debug("\ttotal_paths=%d" % (numpaths))
    if numpaths <= max_exhaustive:
        logging.debug('\tExhaustively enumerating paths')
        for res in find_paths_dfs(G, source, sink, max_exhaustive):
            yield res
    else:
        # return the best path
        logging.debug("\tReturning best path")
        dynprog_reverse(G, sink)
        yield traceback(G, source, sink)
        # return randomly sampled paths
        logging.debug('\tRandomly sampling paths')
        for res in find_paths_markov(G, source, sink, max_exhaustive-1):
            yield res
        # remove graph attributes
        clear_path_attributes(G) 
