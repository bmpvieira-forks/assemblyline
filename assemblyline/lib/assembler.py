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

EDGE_WEIGHT = 'weight'
EDGE_DENSITY = 'density'
EDGE_OUT_FRAC = 'outfrac'
EDGE_IN_FRAC = 'infrac'

PLEN = 'plen'
PWEIGHT = 'pweight'
PDENSITY = 'pdensity'
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
        # add up coverage of edges and set this to the 'score'
        attr_dict = {EDGE_WEIGHT: sum(x.score for x in d['data'])}
        GG[strand].add_edge(u, v, attr_dict=attr_dict)   
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

def visit_node_from_parent(G, parent, child):    
    # compute length, coverage, and rpkm of the current
    # path when extended to include the child
    pattrs = G.node[parent]
    # new path length is parent length + child length    
    path_length = pattrs.get(PLEN, 0) + (child.end - child.start)    
    # to compute the coverage going from the parent to the child, take
    # the existing path coverage at the parent node and multiply it by
    # the outgoing fraction, then add the child coverage multiplied by
    # its incoming fraction.
    parent_cov = pattrs.get(PWEIGHT, pattrs[NODE_WEIGHT])
    child_cov = G.node[child].get(NODE_WEIGHT, 0)
    edata = G.edge[parent][child]    
    path_cov = ((parent_cov * edata[EDGE_OUT_FRAC]) + 
                (child_cov * edata[EDGE_IN_FRAC]))
    # score is total coverage divided by the length of the path
    path_density = path_cov / path_length
    return path_length, path_cov, path_density

def dyn_prog_search(G, source):
    """Find the highest scoring path by dynamic programming"""
    # Taken from NetworkX source code http://networkx.lanl.gov
    stack = [source]
    while stack:
        parent = stack[0]
        for child in G.successors_iter(parent):
            logging.debug('Visiting %s -> %s' % (parent, child))            
            path_length, path_cov, path_density = visit_node_from_parent(G, parent, child)
            logging.debug('\tlen=%d cov=%f density=%f child_density=%f' %
                          (path_length, path_cov, path_density, G.node[child].get(PDENSITY,-1)))
            # only update this node if it has not been visited
            # or if the new score is higher than the old score
            cattrs = G.node[child]
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

def find_best_path(G, source, sink):
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

def recalculate_edge_attrs(G, u, v):
    # find total weight leaving node 'u'
    u_succs = G.successors(u)
    out_weight_total = sum(G.edge[u][x][EDGE_WEIGHT] for x in u_succs)
    if out_weight_total == 0:
        out_weight_total = 1
    # compute edge attributes leaving node 'u'
    for succ in u_succs:
        eattrs = G.edge[u][succ]
        out_frac = eattrs[EDGE_WEIGHT]/float(out_weight_total)
        if out_frac == 0:
            density = 0
        else:
            density = (G.node[u][NODE_WEIGHT] * out_frac)/float(u.end - u.start)
        eattrs[EDGE_OUT_FRAC] = out_frac
        eattrs[EDGE_DENSITY] = density
    # find total density entering node 'v'
    v_preds = G.predecessors(v)
    in_density_total = sum(G.edge[x][v][EDGE_DENSITY] for x in v_preds)
    if in_density_total == 0:
        in_density_total = 1
    # compute edge density fractions entering node 'v'
    for pred in v_preds:
        eattrs = G.edge[pred][v]
        in_frac = eattrs[EDGE_DENSITY]/float(in_density_total)
        eattrs[EDGE_IN_FRAC] = in_frac

def subtract_path(G, path, path_density):
    '''
    subtract the coverage density incurred by traversing the 
    given path from each node and edge in the graph
    '''
    for i,v in enumerate(path):
        # subtract coverage from node
        v_cov = path_density * (v.end - v.start)
        G.node[v][NODE_WEIGHT] = max(0, G.node[v][NODE_WEIGHT] - v_cov)
        if i == 0:
            continue
        # subtract coverage from edge
        u = path[i-1]
        u_cov = path_density * (u.end - u.start)
        eattrs = G.edge[u][v]
        eattrs[EDGE_WEIGHT] = max(0, eattrs[EDGE_WEIGHT] - u_cov)
        # recompute edge attrs
        recalculate_edge_attrs(G, u, v)

def add_dummy_start_end_nodes(G, start_nodes, end_nodes):
    # add 'dummy' tss nodes if necessary
    G.add_node(DUMMY_START_NODE, weight=0)
    for start_node in start_nodes:        
        logging.debug('adding dummy %s -> %s' % (DUMMY_START_NODE, start_node))
        attr_dict = {EDGE_WEIGHT: 0.0,
                     EDGE_DENSITY: 0.0,
                     EDGE_OUT_FRAC: 0.0,
                     EDGE_IN_FRAC: 1.0}                     
        G.add_edge(DUMMY_START_NODE, start_node, attr_dict=attr_dict)
    # add a single 'dummy' end node
    G.add_node(DUMMY_END_NODE, weight=0)
    for end_node in end_nodes:
        logging.debug('adding dummy %s -> %s' % (end_node, DUMMY_END_NODE))
        attr_dict = {EDGE_WEIGHT: 0.0,
                     EDGE_DENSITY: 0.0,
                     EDGE_OUT_FRAC: 1.0,
                     EDGE_IN_FRAC: 0.0}
        G.add_edge(end_node, DUMMY_END_NODE, attr_dict=attr_dict)

def calculate_edge_attrs(G):
    for u,nbrdict in G.adjacency_iter():
        # find total coverage flowing out of this node
        out_score_total = sum(eattrs[EDGE_WEIGHT] for eattrs in nbrdict.itervalues())        
        for v, eattrs in nbrdict.iteritems():
            # find fraction flowing out of each node
            out_frac = eattrs[EDGE_WEIGHT]/float(out_score_total)
            # compute coverage density flowing out of edge
            weight = (G.node[u][NODE_WEIGHT] * out_frac)/float(u.end - u.start)
            eattrs[EDGE_OUT_FRAC] = out_frac
            eattrs[EDGE_DENSITY] = weight
    # now fraction of weight flowing in to each node
    # reverse edge directions and use adjacency iter again
    G.reverse(copy=False)
    for u,nbrdict in G.adjacency_iter():
        # find total coverage density flowing into this node
        in_weight_total = sum(eattrs[EDGE_DENSITY] for eattrs in nbrdict.itervalues())
        for v, eattrs in nbrdict.iteritems():
            # find fraction flowing into node            
            in_frac = eattrs[EDGE_DENSITY]/float(in_weight_total)
            eattrs[EDGE_IN_FRAC] = in_frac
    # reverse the edges back to normal
    G.reverse(copy=False)

def find_suboptimal_paths(G, start_nodes, end_nodes, 
                          fraction_major_path,
                          max_paths):
    # compute initial edge weights and coverage flow in/out of nodes
    calculate_edge_attrs(G)    
    # at an artificial node at the start and end so that all start
    # nodes are searched together for the best path
    add_dummy_start_end_nodes(G, start_nodes, end_nodes)
    score_limit = -1
    while True:        
        # find highest scoring path through graph
        path, path_length, path_cov, path_density = \
            find_best_path(G, DUMMY_START_NODE, DUMMY_END_NODE)
        logging.debug("Path density=%f limit=%f cov=%.2f length=%d path=%s" % 
                      (path_density, score_limit, path_cov, path_length, path))
        # set score limit if this is the first path
        if score_limit < 0:
            score_limit = path_density * fraction_major_path
        elif path_density <= score_limit:
            break
        yield path_density, path[1:-1]
        # subtract the coverage of this path from the graph and recompute 
        # the edge weights
        subtract_path(G, path[1:-1], path_density)
        # remove path node attributes from graph before
        # calling dynamic programming algorithm again
        clear_path_attributes(G)
    # remove dummy nodes from graph
    G.remove_node(DUMMY_START_NODE)
    G.remove_node(DUMMY_END_NODE)

def assemble_subgraph(G, strand, fraction_major_path, max_paths):
    # find start and end nodes in graph
    tss_id = 0
    start_nodes, end_nodes = find_start_and_end_nodes(G, strand)
    logging.debug("START NODES (TSSs): %s" % start_nodes)
    logging.debug("END NODES: %s" % end_nodes)
    for density, path in find_suboptimal_paths(G, start_nodes, end_nodes,
                                               fraction_major_path=fraction_major_path,
                                               max_paths=max_paths):
        # get tss_id from path
        tss_id = G.node[path[0]][NODE_TSS_ID]
        yield tss_id, 1.0e3 * density, path


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
