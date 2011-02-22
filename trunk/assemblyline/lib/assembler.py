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
NODE_WEIGHT = 'weight'
NODE_TSS_ID = 'tss_id'

EDGE_SCORE = 'score'
#EDGE_WEIGHT = 'weight'
EDGE_OUT_FRAC = 'outfrac'
EDGE_IN_FRAC = 'infrac'

PLEN = 'plen'
PCOV = 'pcov'
PFRAC = 'pfrac'
PSCORE = 'pscore'
PSRC = 'psrc'

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

def subtract_path_edge(u, v, nbrdict, path_frac):
    # find total coverage flowing out of this node
    total_weight = sum(eattrs[EDGE_WEIGHT] for eattrs in nbrdict.itervalues())
    # compute coverage used
    used_weight = path_frac * total_weight
    # subtract from path edge
    eattrs = nbrdict[v]    
    eattrs[EDGE_WEIGHT] = max(0, eattrs[EDGE_WEIGHT] - used_weight)
    # subtract from total weight
    total_weight = max(0, total_weight - used_weight)
    # find fraction flowing out of each node
    for v, eattrs in nbrdict.iteritems():
        if total_weight == 0:
            eattrs[EDGE_FRAC] = 0
        else:
            eattrs[EDGE_FRAC] = eattrs[EDGE_WEIGHT]/float(total_weight)

def calculate_edge_weights(G):
    for u,nbrdict in G.adjacency_iter():
        # find total coverage flowing out of this node
        total_weight = sum(eattrs[EDGE_WEIGHT] for eattrs in nbrdict.itervalues())
        # find fraction flowing out of each node
        for v, eattrs in nbrdict.iteritems():
            eattrs[EDGE_FRAC] = eattrs[EDGE_WEIGHT]/float(total_weight)

def calculate_edge_weights(G):
    for u,nbrdict in G.adjacency_iter():
        # find total coverage flowing out of this node
        out_score_total = sum(eattrs[EDGE_SCORE] for eattrs in nbrdict.itervalues())
        for v, eattrs in nbrdict.iteritems():
            # find fraction flowing out of each node
            out_frac = eattrs[EDGE_SCORE]/float(out_score_total)
            # compute coverage density flowing out of edge
            weight = (G.node[u][NODE_WEIGHT] * out_frac)/float(u.end - u.start)
            eattrs[EDGE_OUT_FRAC] = out_frac
            eattrs[EDGE_WEIGHT] = weight
    # now fraction of weight flowing in to each node
    # reverse edge directions and use adjacency iter again
    G.reverse(copy=False)
    for u,nbrdict in G.adjacency_iter():
        # find total coverage density flowing into this node
        in_weight_total = sum(eattrs[EDGE_WEIGHT] for eattrs in nbrdict.itervalues())
        for v, eattrs in nbrdict.iteritems():
            # find fraction flowing into node
            in_frac = eattrs[EDGE_WEIGHT]/float(in_weight_total)
            eattrs[EDGE_IN_FRAC] = in_frac
    
    

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
        # determine strand based on coordinates of nodes
        strand = int(u.start >= v.end)
        # check against strand attribute
        strandattr = d['strand']
        assert cmp_strand(strand, strandattr)
        # add up coverage of edges and set this to the 'weight'
        edge_data_list = d['data']
        # add edge
        attr_dict = {EDGE_SCORE: sum(x.score for x in edge_data_list)}
        GG[strand].add_edge(u, v, attr_dict=attr_dict) 
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
    # annotate each TSS with a unique id
    tss_id = 0
    start_nodes = []
    for tss_nodes in tss_node_dict.itervalues():
        for n in tss_nodes:
            G.node[n][NODE_TSS_ID] = tss_id
            start_nodes.append(n)
        tss_id += 1
    end_nodes = [n for (n,d) in G.out_degree_iter()
                 if (d == 0)]
    return start_nodes, end_nodes

def add_dummy_start_end_nodes(G, start_nodes, end_nodes):
    # add 'dummy' tss nodes if necessary
    total_weight = sum(G.node[n][NODE_WEIGHT] for n in start_nodes)    
    G.add_node(DUMMY_START_NODE, weight=0)
    for start_node in start_nodes:        
        weight = G.node[start_node][NODE_WEIGHT]
        attr_dict = {EDGE_WEIGHT: weight, EDGE_FRAC: (weight / total_weight)}
        logging.debug('adding dummy %s -> %s weight=%f total_weight=%f' % (DUMMY_START_NODE, start_node, weight, total_weight))
        G.add_edge(DUMMY_START_NODE, start_node, attr_dict=attr_dict)
    # add a single 'dummy' end node
    G.add_node(DUMMY_END_NODE, weight=0)
    for end_node in end_nodes:
        weight = G.node[end_node][NODE_WEIGHT]
        attr_dict = {EDGE_WEIGHT: weight, EDGE_FRAC: 1.0}
        logging.debug('adding dummy %s -> %s weight=%f total_weight=%f' % (end_node, DUMMY_END_NODE, total_weight, total_weight))
        G.add_edge(end_node, DUMMY_END_NODE, attr_dict=attr_dict)

def visit_node_from_parent(G, parent, child):    
    # compute length, coverage, and rpkm of the current
    # path when extended to include the child
    pattrs = G.node[parent]
    # new path length is parent length + child length    
    path_length = pattrs.get(PLEN, 0) + (child.end - child.start)    
    # new path coverage is parent coverage + child coverage
    path_cov = pattrs.get(PCOV, 0) + G.node[child].get(NODE_WEIGHT, 0)
    # keep track of the minimum fraction of coverage that
    # flows between and two nodes along the entire path     
    parent_frac = pattrs.get(PFRAC, 1.0)
    child_frac = G.edge[parent][child][EDGE_FRAC]
    path_frac = parent_frac if parent_frac < child_frac else child_frac
    # score is total coverage multiplied by the path fraction
    # divided by the length of the path
    path_score = (1.0e3 * path_frac * path_cov) / path_length
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

                                                
def subtract_path(G, path_frac, path):
    '''
    subtract a fraction of the coverage in each node from the node
    weights and edge scores, then recompute the edge weights
    '''
    for i,v in enumerate(path):
        # recompute node weight
        nattrs = G.node[v]
        used_weight = path_frac * nattrs[NODE_WEIGHT]
        nattrs[NODE_WEIGHT] = max(0, nattrs[NODE_WEIGHT] - used_weight)    
        if i == 0:
            continue
        # recompute edge weight
        u = path[i-1]
        subtract_path_edge(u, v, G.edge[u], path_frac)


def find_suboptimal_paths(G, start_node, end_node, 
                          fraction_major_path,
                          max_paths):
    # make a copy of the graph
    # G = H.copy()
    # iteratively produce paths
    score_limit = -1
    while True:        
        # find highest scoring path through graph
        path, path_score, path_frac, path_cov, path_length = \
            find_best_path(G, start_node, end_node)
        logging.debug("Path score=%f limit=%f frac=%f cov=%.2f length=%d path=%s" % 
                      (path_score, score_limit, path_frac, path_cov, path_length, path))
        # set score limit if this is the first path
        if score_limit < 0:
            score_limit = path_score * fraction_major_path
        elif path_score < score_limit:
            break
        yield path_score, path
        # subtract the coverage of this path from the graph and recompute 
        # the edge weights
        subtract_path(G, path_frac, path)
        for u,v,d in G.edges_iter(data=True):
            logging.debug("Edge (%s,%s) weight=%f frac=%f" % (u,v,d[EDGE_WEIGHT],d[EDGE_FRAC]))
        

def assemble_subgraph(G, strand, fraction_major_path, max_paths):
    # find start and end nodes in graph
    tss_id = 0
    start_nodes, end_nodes = find_start_and_end_nodes(G, strand)
    logging.debug("START NODES (TSSs): %s" % start_nodes)
    logging.debug("END NODES: %s" % end_nodes)
    # at an artificial node at the start and end so that all start
    # nodes are searched together for the best path
    add_dummy_start_end_nodes(G, start_nodes, end_nodes)
    for score, path in find_suboptimal_paths(G, DUMMY_START_NODE, DUMMY_END_NODE,
                                             fraction_major_path=fraction_major_path,
                                             max_paths=max_paths):
        # get tss_id from path
        tss_id = G.node[path[1]][NODE_TSS_ID]
        yield tss_id, score, path[1:-1]
    # remove dummy nodes from graph
    G.remove_node(DUMMY_START_NODE)
    G.remove_node(DUMMY_END_NODE)

def assemble_transcript_graph(G, fraction_major_path, max_paths):
    if fraction_major_path <= 0:
        fraction_major_path = 1e-8
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
                                                         fraction_major_path=fraction_major_path, 
                                                         max_paths=max_paths):
                yield strand, gene_id, tss_id, score, path
            gene_id += 1
