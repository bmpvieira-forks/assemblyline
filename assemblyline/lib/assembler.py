'''
Created on Feb 20, 2011

@author: mkiyer
'''
import collections
import logging
import networkx as nx

from base import Exon, POS_STRAND, NEG_STRAND, NO_STRAND, cmp_strand

# constants for fake 'start' and 'end' nodes
DUMMY_START_NODE = Exon(-1,-1)
DUMMY_END_NODE = Exon(-2,-2)

# constants for use as graph attributes
PLEN = 'length'
PCOV = 'weight'
PFRAC = 'frac'
PSCORE = 'score'
PSRC = 'source'
EDGE_WEIGHT = 'weight'

def sum_node_weights(exon_data_list):
    '''
    sum exon data scores by strand
    '''
    # sum coverage of node by strand
    strand_scores = [0, 0]
    for edata in exon_data_list:
        strand_scores[POS_STRAND] += edata.scores[POS_STRAND]
        strand_scores[NEG_STRAND] += edata.scores[NEG_STRAND]
    return strand_scores

def calculate_edge_weights(G):
    for u,nbrdict in G.adjacency_iter():
        # find total coverage flowing out of this node
        total_score = 0
        for v, eattrs in nbrdict.iteritems():
            total_score += eattrs['score']
        for v, eattrs in nbrdict.iteritems():
            eattrs['weight'] = eattrs['score']/float(total_score)

def build_strand_specific_graphs(G):
    '''
    build separate subgraphs of G - a forward strand graph and a 
    reverse strand graph (if necessary)
    
    nodes in the graph will have a 'weight' attribute equal to the
    total strand-specific coverage for that node

    edges (u,v) in the graph will have a 'weight' attribute corresponding
    to the fraction of total coverage flowing out of node 'u' that goes 
    to 'v'
    '''
    # make graphs
    GG = (nx.DiGraph(), nx.DiGraph())
    # add nodes first
    for n,ndict in G.nodes_iter(data=True):
        strand_weights = sum_node_weights(ndict['data'])
        for strand,w in enumerate(strand_weights):
            if w > 0:
                GG[strand].add_node(n, weight=w)
    # add edges
    for u,v,d in G.edges_iter(data=True):
        # get strand attribute
        strand1 = d['strand']
        # determine strand based on coordinates of nodes
        strand2 = int(u.start >= v.end)
        assert strand1 == strand2
        # add up coverage of edges and call this the 'score'
        edge_data_list = d['data']
        score = sum(x.score for x in edge_data_list)            
        # add edge
        GG[strand].add_edge(u, v, score=score)
    calculate_edge_weights(GG[POS_STRAND])
    calculate_edge_weights(GG[NEG_STRAND])    
    return GG

def find_start_and_end_nodes(G, strand):
    # find unique starting positions and their 
    # corresponding nodes in the graph
    tss_node_dict = collections.defaultdict(lambda: [])
    for n,d in G.in_degree_iter():
        if d == 0:
            tss_pos = n.end if strand == NEG_STRAND else n.start
            tss_node_dict[tss_pos].append(n)
            print 'START', 'STRAND', strand, 'NODE', n 
    end_nodes = [n for (n,d) in G.out_degree_iter()
                 if (d == 0)]
    return tss_node_dict.values(), end_nodes

def add_dummy_start_end_nodes(G, start_nodes, end_nodes):
    # add 'dummy' tss nodes if necessary
    total_weight = sum(G.node[n]['weight'] for n in start_nodes)    
    G.add_node(DUMMY_START_NODE)
    for start_node in start_nodes:        
        weight = G.node[start_node]['weight']
        logging.debug('adding dummy %s -> %s weight=%f total_weight=%f' % (DUMMY_START_NODE, start_node, weight, total_weight))
        G.add_edge(DUMMY_START_NODE, start_node, weight=(weight/total_weight))
    # add a single 'dummy' end node
    G.add_node(DUMMY_END_NODE)
    for end_node in end_nodes:
        logging.debug('adding dummy %s -> %s weight=%f total_weight=%f' % (end_node, DUMMY_END_NODE, total_weight, total_weight))
        G.add_edge(end_node, DUMMY_END_NODE, weight=1.0)

def visit_node_from_parent(G, parent, child):    
    # compute length, coverage, and rpkm of the current
    # path when extended to include the child
    pattrs = G.node[parent]
    # new path length is parent length + child length    
    path_length = pattrs.get(PLEN, 0) + (child.end - child.start)    
    # new path coverage is parent coverage + child coverage
    path_cov = pattrs.get(PCOV, 0) + G.node[child].get(PCOV, 0)
    # keep track of the minimum fraction of coverage that
    # flows between and two nodes along the entire path     
    parent_frac = pattrs.get(PFRAC, 1.0)
    child_frac = G.edge[parent][child][EDGE_WEIGHT]
    path_frac = parent_frac if parent_frac < child_frac else child_frac
    # score is total coverage multiplied by the path fraction
    # divided by the length of the path
    path_score = (path_frac * path_cov) / path_length
    return path_length, path_cov, path_frac, path_score

def dyn_prog_search(G, source):
    """Find the highest scoring path by dynamic programming"""
    # Taken from NetworkX source code http://networkx.lanl.gov
    stack = [source]
    while stack:
        parent = stack[0]
        for child in G.successors_iter(parent):
            logging.debug('Visiting %s -> %s' % (parent, child))            
            path_length, path_cov, path_frac, path_score = visit_node_from_parent(G, parent, child)
            logging.debug('\tlen=%d frac=%f score=%f child_score=%f' %
                          (path_length, path_frac, path_score, G.node[child].get(PSCORE,-1)))
            # only update this node if it has not been visited
            # or if the new score is higher than the old score
            cattrs = G.node[child]
            if path_score > cattrs.get(PSCORE, -1):
                # keep pointer to parent node that produced this high scoring path
                cattrs[PLEN] = path_length
                cattrs[PCOV] = path_cov
                cattrs[PFRAC] = path_frac
                cattrs[PSCORE] = path_score
                cattrs[PSRC] = parent
                logging.debug('\tupdated child len=%d cov=%f frac=%f score=%f' % 
                              (path_length, path_cov, path_frac, path_score))
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
    frac = pattrs[PFRAC]    
    length = pattrs[PLEN]
    cov = pattrs[PCOV]
    while path[-1] != source:
        #logging.debug('path=%s node=%s attrs=%s' %
        #              (path, path[-1], G.node[path[-1]]))
        path.append(G.node[path[-1]][PSRC])
    path.reverse()
    #logging.debug("FINAL path=%s" % (path))
    # clear path attributes
    for n in G.nodes_iter():
        nattrs = G.node[n]
        if PLEN in nattrs:
            del nattrs[PLEN]
            del nattrs[PCOV]
            del nattrs[PFRAC]
            del nattrs[PSCORE]
            del nattrs[PSRC]
    return path, score, frac, cov, length


def find_suboptimal_paths(G, start_node, end_node, 
                          fraction_major_path,
                          max_paths,
                          max_iters=10000):
    # copy the EDGEWT attribute rather than the entire graph
    edge_weight_dict = {}
    for u,v,d in G.edges_iter(data=True):
        if EDGE_WEIGHT in d:
            edge_weight_dict[(u,v)] = d[EDGE_WEIGHT]
    # find overall highest scoring path through graph
    path, path_score, path_frac, path_cov, path_length = \
        find_best_path(G, start_node, end_node)
    covscore = 1.0e3 * path_cov * path_frac / path_length
    score_limit = covscore * fraction_major_path
    logging.debug("Maximum path covscore=%f limit=%f weight=%f score=%f path=%s" % 
                  (covscore, score_limit, path_frac, path_score, path))
    # enumerate paths until the score falls below the 
    # specified percentage of the best path
    while covscore > score_limit:
        yield covscore, path
        # subtract path weight from edge weights
        for parent,child in zip(path[:-1], path[1:]):
            pcedge = G.edge[parent][child]
            if EDGE_WEIGHT in pcedge:
                pcedge[EDGE_WEIGHT] -= path_frac
        # get next path        
        path, path_score, path_frac, path_cov, path_length = \
            find_best_path(G, start_node, end_node)
        covscore = 1.0e3 * path_cov * path_frac / path_length
        logging.debug("Current path covscore=%f limit=%f weight=%f score=%f path=%s" % 
                      (covscore, score_limit, path_frac, path_score, path))
    # restore original edge weights
    for u,v,d in G.edges_iter(data=True):
        if EDGE_WEIGHT in d:
            d[EDGE_WEIGHT] = edge_weight_dict[(u,v)]

def assemble_subgraph(G, strand, fraction_major_isoform, max_paths, max_iters):
    # find start and end nodes in graph
    tss_id = 0
    start_node_lists, end_nodes = find_start_and_end_nodes(G, strand)
    for start_nodes in start_node_lists:
        logging.debug("START NODES (TSSs): %s" % start_nodes)
        logging.debug("END NODES: %s" % end_nodes)
        # at an artificial node at the start and end that all paths with
        # the same tss are searched together
        add_dummy_start_end_nodes(G, start_nodes, end_nodes)
        for score, path in find_suboptimal_paths(G, DUMMY_START_NODE, DUMMY_END_NODE,
                                                 fraction_major_path=fraction_major_isoform,
                                                 max_paths=max_paths,
                                                 max_iters=max_iters):
            yield tss_id, score, path[1:-1]
        # remove dummy nodes from graph
        G.remove_node(DUMMY_START_NODE)
        G.remove_node(DUMMY_END_NODE)
        tss_id += 1

def assemble_transcript_graph(G, fraction_major_isoform, max_paths, max_iters):
    # transform transcript graph into strand-specific graphs
    # for forward and reverse strands, and calculate node and
    # edge weights
    GG = build_strand_specific_graphs(G)
    # assemble strands one at a time
    for strand, Gstrand in enumerate(GG):
        logging.debug("STRAND: %s" % strand)
        gene_id = 0        
        # get connected components - unconnected components are considered
        # different genes
        for Gsub in nx.weakly_connected_component_subgraphs(Gstrand):
            for tss_id, score, path in assemble_subgraph(Gsub, strand, 
                                                         fraction_major_isoform=fraction_major_isoform, 
                                                         max_paths=max_paths,
                                                         max_iters=max_iters):
                yield strand, gene_id, tss_id, score, path
            gene_id += 1
