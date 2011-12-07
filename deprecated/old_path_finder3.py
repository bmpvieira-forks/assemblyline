'''
Created on Dec 4, 2010

@author: mkiyer
'''
import logging
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
PCOV = 'pcov'
PWT = 'pwt'
PSCORE = 'pscore'
PSRC = 'psrc'

EDGEWT = 'psi'
EDGECOV = 'weight'

def visit_node_from_parent(G, parent, child):    
    pattrs = G.node[parent]
    # new path length is parent length + child length    
    parent_path_length = pattrs.get(PLEN, 0)
    path_length = parent_path_length + (child.end - child.start)    
    # new path coverage is parent cov + child cov
    parent_path_cov = pattrs.get(PCOV, 0)
    path_cov = parent_path_cov + G.node[child].get(EDGECOV, 0)
    # new path weight is minimum of parent weight (if exists) and child weight
    pcedge = G.edge[parent][child]
    if EDGEWT in pcedge:
        child_weight = pcedge[EDGEWT]
        parent_weight = pattrs.get(PWT, child_weight)  
    else:
        # if child does not have an EDGEWT field use the
        # parent weight
        parent_weight = pattrs[PWT]   
        child_weight = parent_weight
    # update weight value
    path_weight = child_weight if child_weight < parent_weight else parent_weight
    # score if weight / length
    path_score = path_weight / path_length
    return path_length, path_cov, path_weight, path_score

def dyn_prog_search(G,source):
    """Find the highest scoring path by dynamic programming"""
    # Taken from NetworkX source code http://networkx.lanl.gov
    stack = [source]
    while stack:
        parent = stack[0]
        for child in G.successors_iter(parent):
            logging.debug('Visiting %s -> %s' % (parent, child))
            path_length, path_cov, path_weight, path_score = visit_node_from_parent(G, parent, child)
            logging.debug('\tlen=%d wt=%f score=%f child_score=%f' %
                          (path_length, path_weight, path_score, G.node[child].get(PSCORE,-1)))
            cattrs = G.node[child]
            # only update this node if it has not been visited
            # or if the new score is higher than the old score
            if (PSCORE not in cattrs) or (path_score > cattrs[PSCORE]):
                # keep pointer to parent node that produced this high scoring path
                cattrs[PLEN] = path_length
                cattrs[PCOV] = path_cov
                cattrs[PWT] = path_weight
                cattrs[PSCORE] = path_score
                cattrs[PSRC] = parent
                logging.debug('\tupdated child len=%d cov=%f wt=%f score=%f' % 
                              (path_length, path_cov, path_weight, path_score))
                # continue iterating through the child path
                stack.append(child)
        stack.pop(0)

def find_best_path(G, source, sink):
    '''
    use dynamic programming to find the highest scoring path through 
    the graph starting from 'source'
    '''
    # run the dynamic programming algorithm to
    # score paths
    dyn_prog_search(G, source)
    # traceback    
    path = [sink]
    pattrs = G.node[sink]
    score = pattrs[PSCORE]
    weight = pattrs[PWT]    
    length = pattrs[PLEN]
    cov = pattrs[PCOV]
    while path[-1] != source:
        logging.debug('path=%s node=%s attrs=%s' %
                      (path, path[-1], G.node[path[-1]]))
        path.append(G.node[path[-1]][PSRC])
    path.reverse()
    logging.debug("FINAL path=%s" % (path))
    # clear path attributes
    for n in G.nodes_iter():
        nattrs = G.node[n]
        if PLEN in nattrs:
            del nattrs[PLEN]
            del nattrs[PCOV]
            del nattrs[PWT]
            del nattrs[PSCORE]
            del nattrs[PSRC]
    return path, score, weight, cov, length

def find_suboptimal_paths(G, start_node, end_node, 
                          fraction_major_path,
                          max_paths,
                          max_iters=10000):
    # copy the EDGEWT attribute rather than the entire graph
    edge_weight_dict = {}
    for u,v,d in G.edges_iter(data=True):
        if EDGEWT in d:
            edge_weight_dict[(u,v)] = d[EDGEWT]
    # make a copy of the graph
    #G = H.copy()
    # get edges in breadth-first order for dynamic programming  
    path, path_score, path_weight, path_cov, path_length = find_best_path(G, start_node, end_node)
    covscore = 1.0e3 * path_cov * path_weight / path_length
    score_limit = covscore * fraction_major_path
    logging.debug("Maximum path covscore=%f limit=%f weight=%f psiscore=%f path=%s" % (covscore, score_limit, path_weight, path_score, path))
    # enumerate paths until the score falls below the 
    # specified percentage of the best path
    while covscore > score_limit:
        yield covscore, path
        # subtract path weight from edge weights
        for parent,child in zip(path[:-1], path[1:]):
            pcedge = G.edge[parent][child]
            if EDGEWT in pcedge:
                pcedge[EDGEWT] -= path_weight
        # get next path        
        path, path_score, path_weight, path_cov, path_length = find_best_path(G, start_node, end_node)
        covscore = 1.0e3 * path_cov * path_weight / path_length
        logging.debug("Current path covscore=%f limit=%f weight=%f psiscore=%f path=%s" % (covscore, score_limit, path_weight, path_score, path))
    # restore original edge weights
    for u,v,d in G.edges_iter(data=True):
        if EDGEWT in d:
            d[EDGEWT] = edge_weight_dict[(u,v)]

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
            H.add_node(n, weight=weight, ids=attr_dict['ids'])
    # add exon-exon edges
    for n,attr_dict in G.nodes_iter(data=True):
        if n.node_type == INTRON:
            weight = sum(id_score_map[id] for id in attr_dict['ids'])            
            for pred in G.predecessors_iter(n):
                for succ in G.successors_iter(n):   
                    psi = (weight/H.node[pred]['weight'])               
                    H.add_edge(pred, succ, weight=weight, 
                               psi=psi, ids=attr_dict['ids'])
    return H

def find_start_and_end_nodes(G):
    # find graph strand
    strand = reduce(merge_strand, iter(n.strand for n in G.nodes_iter()))
    # find unique starting positions and their 
    # corresponding nodes in the subgraph
    tss_node_dict = collections.defaultdict(lambda: [])
    for n,d in G.in_degree_iter():
        if d == 0:
            tss_pos = n.end if strand == NEG_STRAND else n.start
            tss_node_dict[tss_pos].append(n)
    # add a single 'dummy' end node
    end_nodes = [n for (n,d) in G.out_degree_iter()
                 if (d == 0)]
    logging.debug("START NODES (TSSs): %s" % tss_node_dict.values())
    logging.debug("END NODES: %s" % end_nodes)
    return tss_node_dict.values(), end_nodes

def add_dummy_start_end_nodes(G, start_nodes, end_nodes):
    # add 'dummy' tss nodes if necessary
    dummy_start_node = Node(0, 0, NO_STRAND, DUMMY)
    G.add_node(dummy_start_node)
    total_weight = sum(G.node[n]['weight'] for n in start_nodes)
    for start_node in start_nodes:
        weight = G.node[start_node]['weight']
        psi = weight / total_weight
        logging.debug('adding dummy %s -> %s psi=%f' % (dummy_start_node, start_node, psi))
        G.add_edge(dummy_start_node, start_node, 
                   weight=G.node[start_node]['weight'], 
                   psi=psi, ids=[])            
    # add a single 'dummy' end node
    dummy_end_node = Node(-1, -1, NO_STRAND, DUMMY)
    G.add_node(dummy_end_node)
    for end_node in end_nodes:
        G.add_edge(end_node, dummy_end_node)    
    return dummy_start_node, dummy_end_node

def get_isoforms(G, transcripts,
                 fraction_major_path=0.15,
                 max_paths=10, 
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
        start_node_lists, end_nodes = find_start_and_end_nodes(Hsubgraph)
        for start_nodes in start_node_lists:
            logging.debug('start nodes %s' % (start_nodes))
            dummy_start_node, dummy_end_node = add_dummy_start_end_nodes(Hsubgraph, start_nodes, end_nodes)            
            for score, path in find_suboptimal_paths(Hsubgraph, dummy_start_node, dummy_end_node, 
                                                     fraction_major_path=fraction_major_path,
                                                     max_paths=max_paths, 
                                                     max_iters=max_iters):
                # get full set of transcripts that support this path
                
                # remove dummy nodes (first and last) when returning path
                yield gene_id, tss_id, score, path[1:-1]
            # remove dummy nodes from graph
            Hsubgraph.remove_node(dummy_start_node)
            Hsubgraph.remove_node(dummy_end_node)
            tss_id += 1
        gene_id += 10
    logging.debug("/PATHFINDER")
