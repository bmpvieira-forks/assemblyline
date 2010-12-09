'''
Created on Dec 4, 2010

@author: mkiyer
'''
import logging
import networkx as nx
import collections
import heapq
import operator

from bx.intervals.intersection import Interval, IntervalTree
from bx.intervals.cluster import ClusterTree

from cNode import cmp_strand, strand_str_to_int
from isoform_graph import POS_STRAND, NEG_STRAND, NO_STRAND, INTRON, EXON

def imin2(x,y):
    return x if x <= y else y
def imax2(x,y):
    return x if x >= y else y

def rank_paths(G, parent, length, weight):
    children = G.successors_iter(parent)
    pedge = G.edge[parent]
    # score is the minimum of the current weight and the 
    # weight of the edge, normalized to transcript length
    h = [(-imin2(weight,pedge[c]['weight'])/float(length + (c.end - c.start)), c) for c in children]    
    #h = [(-(weight + pedge[c]['weight'])/float(length + (c.end - c.start)), c) for c in children]
    #h = [(-((weight + G.node[c]['weight'])/float(length + (c.end - c.start))), c) for c in children]
    heapq.heapify(h)
    while len(h) > 0:
        score, child = heapq.heappop(h)
        print 'score', score, 'child', child
        yield child

#def dfs_paths(G, start_node, end_nodes):
#    """
#    Produce edges in a depth-first-search starting at source
#    
#    Based on http://www.ics.uci.edu/~eppstein/PADS/DFS.py
#    by D. Eppstein, July 2004.
#    Modified from the code in networkx http://networkx.lanl.gov
#    """
#    # produce edges for components with source
#    path = [start_node]
#    path_length = start_node.end - start_node.start
#    path_weight = G.node[start_node]['weight']
#    stack = [(path[-1], rank_paths(G, path[-1], path_length, path_weight))]
#    while stack:
#        parent, children = stack[-1]
#        #print 'current path', str(path), path_length, path_weight
#        try:
#            child = next(children)
#            # add to path
#            path.append(child)
#            path_length += (child.end - child.start)
#            path_weight += G.edge[parent][child]['weight']
#            #path_weight += G.node[child]['weight']
#            # check for finished path
#            if path[-1] in end_nodes:
#                yield list(path), path_length, path_weight           
#            stack.append((path[-1], rank_paths(G, path[-1], path_length, path_weight))) 
#        except StopIteration:
#            print 'popping stack', stack[-1]
#            lastchild,childiter = stack[-1]
#            print 'popping stack', stack[-1]
#            print 'current parent', parent
#            print 'last child', lastchild
#            lastchild = path.pop()
#            print 'last path child', lastchild
#            path_length -= (lastchild.end - lastchild.start)
#            #path_weight -= G.node[lastchild]['weight']
#            stack.pop()

def dfs_paths(G, start_node, end_nodes):
    """
    Produce edges in a depth-first-search starting at source
    
    Based on http://www.ics.uci.edu/~eppstein/PADS/DFS.py
    by D. Eppstein, July 2004.
    Modified from the code in networkx http://networkx.lanl.gov
    """
    # produce edges for components with source
    path_length = start_node.end - start_node.start
    path_weight = G.node[start_node]['weight']
    path = [start_node]
    stack = [(path_weight, rank_paths(G, path[-1], path_length, path_weight))]
    while stack:
        parent = path[-1]
        parent_weight, children = stack[-1]
        print 'current w=%f l=%d p=%s' % (path_weight, path_length, path)
        try:
            child = next(children)
            # add to path
            path_length += (child.end - child.start)
            path_weight = imin2(path_weight, G.edge[parent][child]['weight'])
            path.append(child)
            # check for finished path
            if path[-1] in end_nodes:
                yield list(path), path_length, path_weight           
            stack.append((path_weight, rank_paths(G, path[-1], path_length, path_weight))) 
        except StopIteration:
            stack.pop()
            path.pop()
            path_length -= (parent.end - parent.start)
            path_weight = parent_weight
            print 'revert to w=%f l=%d p=%s' % (path_weight, path_length, path)
            

def find_suboptimal_paths(G, start_node, end_nodes, max_paths=5, max_iters=10000):
    debug_every = 1000
    debug_next = 1000
    iters = 0
    paths = []
    for path, path_length, path_weight in dfs_paths(G, start_node, end_nodes):        
        score = 1.0e3 * path_weight / path_length
        logging.debug("score=%f path=%s" % (score, path))
        if len(paths) < max_paths:
            # add path to graph
            heapq.heappush(paths, (score, path))
        else:
            # add path and discard the lowest scoring path
            heapq.heappushpop(paths, (score, path))
        iters += 1
        if iters == debug_next:
            logging.debug("Path finder examined %d paths" % (iters))
            debug_next += debug_every
        if iters == max_iters:
            logging.warning("Path scoring reached max iterations before enumerating all paths")
            break
    while len(paths) > 0:
        yield heapq.heappop(paths)

def get_transcript_score_map(transcripts):
    # get transcript scores
    # TODO: could do this when originally reading in data
    id_score_map = {}
    for t in transcripts:
        # add scores to the score lookup table
        id_score_map[t.id] = t.score
        for e in t.exons:
            id_score_map[e.id] = e.score
    return id_score_map

def transform_graph(G, id_score_map):
    '''
    convert introns to weighted edges rather than nodes and
    use 'ids' attribute of nodes to retrieve scores for that
    node as well as overall node/edge weights    
    '''
    # create a new digraph
    H = nx.DiGraph()
    # add all nodes from previous graph
    for n,attr_dict in G.nodes_iter(data=True):
        # TODO: keep track of sample-specific scores
        if n.node_type == EXON:
            weight = sum(id_score_map[id] for id in attr_dict['ids'])
            #weight = 1.0e3 * count / float(n.end - n.start)
            H.add_node(n, weight=weight, ids=attr_dict['ids'])
    # add exon-exon edges
    for n,attr_dict in G.nodes_iter(data=True):
        if n.node_type == INTRON:
            weight = sum(id_score_map[id] for id in attr_dict['ids'])
            for pred in G.predecessors_iter(n):
                for succ in G.successors_iter(n):                    
                    H.add_edge(pred, succ, weight=weight, ids=attr_dict['ids'])
    return H

def get_isoforms(G, transcripts):
    # map transcript ids to scores
    id_score_map = get_transcript_score_map(transcripts)
    # replace intron nodes with exon-exon edges
    H = transform_graph(G, id_score_map)
    # partition graph into connected components
    gene_id = 0
    for Hsubgraph in nx.weakly_connected_component_subgraphs(H):
        # find all source nodes in subgraph
        start_nodes = [n for (n,d) in Hsubgraph.in_degree_iter()
                       if (d == 0)]
        end_nodes = [n for (n,d) in Hsubgraph.out_degree_iter()
                     if (d == 0)]
        tss_id = 0
        for start_node in start_nodes:
            for score, path in find_suboptimal_paths(Hsubgraph, start_node, end_nodes, max_paths=5, max_iters=10000):
                yield gene_id, tss_id, score, path
            tss_id += 1
        gene_id += 1
