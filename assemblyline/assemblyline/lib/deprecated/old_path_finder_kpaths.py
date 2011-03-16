'''
Created on Dec 4, 2010

@author: mkiyer
'''
import logging
import collections
import bisect
import networkx as nx

from base import NODE_WEIGHT, EDGE_OUT_FRAC, EDGE_IN_FRAC, EDGE_DENSITY

NodePtr = collections.namedtuple("NodePtr", ["density", "weight", "length", "src"])

PATH_PTRS = 'pptrs'
WORST_DENSITY = 'pworst'

PLEN = 'plen'
PWEIGHT = 'pweight'
PDENSITY = 'pdensity'
PSRC = 'psrc'

def calc_score(G, parent, pattrs, child):
    """`parent` and `child` are :class:`NodePtr` objects
    """
    # new path length is parent length + child length
    length = pattrs.get(PLEN, 0) + (child.end - child.start)
    # to compute the coverage going from the parent to the child, take
    # the existing path coverage at the parent node and multiply it by
    # the outgoing fraction, then add the child coverage multiplied by
    # its incoming fraction.
    child_weight = G.node[child].get(NODE_WEIGHT, 0)
    edata = G.edge[parent][child]
    parent_weight = pattrs.get(PWEIGHT, pattrs[NODE_WEIGHT])
    weight = (parent_weight * edata[EDGE_IN_FRAC] + 
              child_weight * edata[EDGE_OUT_FRAC])
    # score is total coverage divided by the length of the path
    density = weight / length
    return density, weight, length

def dynprog_forward(G, sink):
    """Compute the highest scoring path from each node to the
    end node"""
    # reverse edge directions and use adjacency iter again
    G.reverse(copy=False)
    # Taken from NetworkX source code http://networkx.lanl.gov
    stack = [sink]
    while stack:
        parent = stack[0]
        for child in G.successors_iter(parent):
            logging.debug('Visiting %s -> %s' % (parent, child))
          



def visit_node_from_parent(G, parent, child, max_paths):
    # compute length, weight, and density of the current
    # path when extended to include the child
    pattrs = G.node[parent]
    pptrs = pattrs[PATH_PTRS]
    cattrs = G.node[child]
    cptrs = cattrs[PATH_PTRS]
    for pptr in pptrs:
        density, weight, length = calc_score(G, parent, pptr, child)
        if len(cptrs) < max_paths:
            # fill up the list of NodePtrs until its length equals 'max_paths'
            bisect.insort_left(cptrs, NodePtr(-density, weight, length, parent))
        else:
            worst_density = -(cptrs[-1].density)
            if density > worst_density:
                # replace the worst item with this one
                cptrs.pop()
                bisect.insort_left(cptrs, NodePtr(-density, weight, length, parent))
            else:
                # given that parent pointers are already sorted from best
                # to worst, we will not find anything else to insert here
                break

      

def dyn_prog_search(G, source, max_paths):
    """Find the highest scoring path by dynamic programming"""
    # Taken from NetworkX source code http://networkx.lanl.gov
    stack = [source]
    while stack:
        parent = stack[0]
        for child in G.successors_iter(parent):
            logging.debug('Visiting %s -> %s' % (parent, child))

            visit_node_from_parent(G, parent, child, max_paths):
            # continue iterating through the child path
            stack.append(child)

            # only update this node if it has not been visited
            # or if the new score is higher than the old score
            cattrs = G.node[child]
            cptrs = cattrs[PPTRS]
            num_cptrs = len(cptrs)
            cworst = cattrs[PWORST]
            for nptr in visit_node_from_parent(G, parent, child):

            path_length, path_cov, path_density = visit_node_from_parent(G, parent, child)
            logging.debug('\tlen=%d cov=%f density=%f child_density=%f' %
                          (path_length, path_cov, path_density, G.node[child].get(PDENSITY,-1)))

            if path_density > cattrs.get(PDENSITY, -1):
                # keep pointer to parent node that produced this high scoring path
                cattrs[PLEN] = path_length
                cattrs[PWEIGHT] = path_cov
                cattrs[PDENSITY] = path_density
                cattrs[PSRC] = parent
                logging.debug('\tupdated child len=%d cov=%f density=%f' % 
                              (path_length, path_cov, path_density))
                # continue iterating through the child path
                stack.append(child)
        stack.pop(0)

def find_best_path(G, source, sink, max_paths=100):
    '''
    use dynamic programming to find the highest scoring path through 
    the graph starting from 'source' and ending at 'sink'
    '''
    # run the dynamic programming algorithm to
    # score paths
    dyn_prog_search(G, source)
    # traceback
    path = [sink]
    pattrs = G.node[sink]
    score = pattrs[PDENSITY]
    length = pattrs[PLEN]
    cov = pattrs[PWEIGHT]
    while path[-1] != source:
        #logging.debug('path=%s node=%s attrs=%s' %
        #              (path, path[-1], G.node[path[-1]]))
        path.append(G.node[path[-1]][PSRC])
    path.reverse()
    #logging.debug("FINAL path=%s" % (path))
    return path, length, cov, score
