'''
Created on Dec 4, 2010

@author: mkiyer
'''
import logging
import heapq
import operator
import collections
import networkx as nx

from base import merge_strand, POS_STRAND, NEG_STRAND, NO_STRAND, INTRON, EXON, DUMMY
from cNode import Node

def imin2(x,y):
    return x if x <= y else y
def imax2(x,y):
    return x if x >= y else y

# constants for use as graph attributes
PLEN = 'plen'
PWT = 'pwt'
PSCORE = 'pscore'
PSRC = 'psrc'

def bfs_edges(G,source):
    """Produce edges in a breadth-first-search starting at source."""
    # Taken from NetworkX source code http://networkx.lanl.gov
    visited=set([source])
    stack = [(source,iter(G[source]))]
    while stack:
        parent,children = stack[0]
        try:
            child = next(children)
            if child not in visited:
                yield parent,child
                visited.add(child)
                stack.append((child,iter(G[child])))
        except StopIteration:
            stack.pop(0)

def visit_node_from_parent(G, parent, child):    
    pnode = G.node[parent]
    cnode = G.node[child]
    # new path length is parent length + child length    
    parent_path_length = pnode.get(PLEN, 0)
    path_length = parent_path_length + (child.end - child.start)    
    # new path weight is minimum of parent weight (if exists) and child weight
    pcedge = G.edge[parent][child]
    if 'weight' in pcedge:
        child_weight = pcedge['weight']
        parent_weight = pnode.get(PWT, child_weight)  
    else:
        # if child does not have a 'weight' field use the
        # parent weight
        parent_weight = pnode[PWT]   
        child_weight = parent_weight
    path_weight = child_weight if child_weight < parent_weight else parent_weight
    # score if weight / length
    path_score = path_weight / path_length
    child_path_score = cnode.get(PSCORE, path_score)
    print 'visiting', parent, '->', child, 'len', path_length, 'wt', path_weight, 'score', path_score, 'child score', child_path_score
    if path_score >= child_path_score:
        # keep pointer to parent node that produced this high scoring path
        cnode = G.node[child]
        cnode[PLEN] = path_length
        cnode[PWT] = path_weight
        cnode[PSCORE] = path_score
        cnode[PSRC] = parent
        print 'wrote child', parent, '->', child, 'len', path_length, 'wt', path_weight, 'score', path_score

def clear_path_attributes(G):
    for n in G.nodes_iter():
        G.node[n] = {}


def find_best_path(G, sorted_edges, source, sink):
    '''
    use dynamic programming to find the highest scoring path through 
    the graph starting from 'source'
    '''
    # topologically sort nodes
    for parent, child in sorted_edges:
        visit_node_from_parent(G, parent, child)
    # traceback    
    path = [sink]
    attrdict = G.node[sink]
    score = attrdict[PSCORE]
    weight = attrdict[PWT]
    while path[-1] != source:
        path.append(G.node[path[-1]][PSRC])
    path.reverse()
    return path, score, weight


def rank_paths(G, parent, length, weight):
    children = G.successors_iter(parent)
    pedge = G.edge[parent]
    # score is the minimum of the current weight and the 
    # weight of the edge, normalized to transcript length
    h = [(-imin2(weight,pedge[c]['weight'])/float(length + (c.end - c.start)), c) for c in children]    
    heapq.heapify(h)
    while len(h) > 0:
        yield heapq.heappop(h)[1]

def dfs_paths(G, start_node, end_nodes):
    """
    performs depth-first search beginning at 'start_node',
    and prioritizes paths based on the 'rank_paths' weighting
    function.  the function is greedy in that first enumerates
    all the path extensions from the current path and chooses
    the one that maximizes the path score    
    """
    # produce edges for components with source
    path_length = start_node.end - start_node.start
    path_weight = G.node[start_node]['weight']
    path = [start_node]
    # check if start node is also end node
    if start_node in end_nodes:
        yield list(path), path_length, path_weight
        return
    stack = [(path_weight, rank_paths(G, path[-1], path_length, path_weight))]
    while stack:
        parent = path[-1]
        parent_weight, children = stack[-1]
        #print 'current w=%f l=%d p=%s' % (path_weight, path_length, path)
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
            #print 'revert to w=%f l=%d p=%s' % (path_weight, path_length, path)

#def find_suboptimal_paths(G, start_node, end_nodes, 
#                          fraction_major_path,
#                          max_paths,
#                          max_iters=10000):
#    debug_every = 1000
#    debug_next = 1000
#    iters = 0
#    paths = []
#    for path, path_length, path_weight in dfs_paths(G, start_node, end_nodes):        
#        score = 1.0e3 * path_weight / path_length
#        #logging.debug("score=%f path=%s" % (score, path))
#        if len(paths) < max_paths:
#            # add path to graph
#            heapq.heappush(paths, (score, path))
#        else:
#            # add path and discard the lowest scoring path
#            heapq.heappushpop(paths, (score, path))
#        iters += 1
#        if iters == debug_next:
#            logging.debug("Path finder examined %d paths" % (iters))
#            debug_next += debug_every
#        if iters == max_iters:
#            logging.warning("Path scoring reached max iterations before enumerating all paths")
#            break
#    # sort from hi -> low score
#    paths = sorted(paths, key=operator.itemgetter(0), reverse=True)[:max_paths]
#    path_score_limit = paths[0][0] * fraction_major_path
#    for score, path in paths:
#        if score < path_score_limit:
#            break
#        yield score, path
#    del paths

def find_suboptimal_paths(G, start_node, end_node, 
                          fraction_major_path,
                          max_paths,
                          max_iters=10000):
    # topologically sort nodes
    sorted_nodes = nx.topological_sort(G)
    # get sorted edges of the sorted nodes
    sorted_edges = []
    for n in sorted_nodes:
        sorted_edges.extend((p,n) for p in G.predecessors(n))
    
    H = G.copy()
    path, path_score, path_weight = find_best_path(H, sorted_edges, start_node, end_node)
    score_limit = path_score * fraction_major_path
    # enumerate paths until the score falls below the 
    # specified percentage of the best path
    while path_score >= score_limit:
        # subtract path weight from edge weights
        for parent,child in zip(path[:-1], path[1:]):
            pcedge = G.edge[parent][child]
            if 'weight' in pcedge:
                #print 'weight', parent, '->', child, pcedge['weight']  
                pcedge['weight'] -= path_weight
                #print 'weight after', pcedge['weight']
        #print 'path', path, path_score, path_weight
        yield path_score, path
        H = G.copy()
        path, path_score, path_weight = find_best_path(H, sorted_edges, start_node, end_node)

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

def get_start_and_end_nodes(G):
    # find graph strand
    strand = reduce(merge_strand, iter(n.strand for n in G.nodes_iter()))
    # find unique starting positions and their 
    # corresponding nodes in the subgraph
    tss_node_dict = collections.defaultdict(lambda: [])
    for n,d in G.in_degree_iter():
        if d == 0:
            tss_pos = n.end if strand == NEG_STRAND else n.start
            tss_node_dict[tss_pos].append(n)
    # add 'dummy' tss nodes if necessary
    dummy_start_nodes = []
    for tss_pos, tss_nodes in tss_node_dict.iteritems():
        dummy_start_node = Node(0, 0, NO_STRAND, DUMMY)
        dummy_start_nodes.append(dummy_start_node)
        #weight = sum(G.node[n]['weight'] for n in tss_nodes)
        #G.add_node(dummy_start_node, weight=weight, ids=[])
        G.add_node(dummy_start_node)
        for tss_node in tss_nodes:
            #G.add_edge(dummy_start_node, tss_node)            
            G.add_edge(dummy_start_node, tss_node, weight=G.node[tss_node]['weight'], ids=[])            
    # add a single 'dummy' end node
    end_nodes = [n for (n,d) in G.out_degree_iter()
                 if (d == 0)]
    dummy_end_node = Node(-1, -1, NO_STRAND, DUMMY)
    G.add_node(dummy_end_node)
    #G.add_node(dummy_end_node, weight=weight, ids=[])
    for end_node in end_nodes:
        G.add_edge(end_node, dummy_end_node)    
        #G.add_edge(end_node, dummy_end_node, weight=0, ids=[])    
    return dummy_start_nodes, dummy_end_node

def get_isoforms(G, transcripts,
                 fraction_major_path=0.15,
                 max_paths=5, 
                 max_iters=10000):                 
    # map transcript ids to scores
    id_score_map = get_transcript_score_map(transcripts)
    # replace intron nodes with exon-exon edges
    H = transform_graph(G, id_score_map)
    # partition graph into connected components
    gene_id = 0
    tss_id = 0
    logging.debug("PATHFINDER")
    for Hsubgraph in nx.weakly_connected_component_subgraphs(H):
        # add dummy start and end nodes to the graph
        dummy_start_nodes, dummy_end_node = get_start_and_end_nodes(Hsubgraph)
        for start_node in dummy_start_nodes:
            for score, path in find_suboptimal_paths(Hsubgraph, start_node, dummy_end_node, 
                                                     fraction_major_path=fraction_major_path,
                                                     max_paths=max_paths, 
                                                     max_iters=max_iters):
                # remove dummy nodes when returning path
                yield gene_id, tss_id, score, path[1:-1]
            tss_id += 1
        gene_id += 1
    logging.debug("/PATHFINDER")
