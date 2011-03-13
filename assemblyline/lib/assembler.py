'''
Created on Feb 20, 2011

@author: mkiyer
'''
import collections
import logging
import networkx as nx

from bx.cluster import ClusterTree

from base import Exon, POS_STRAND, NEG_STRAND, NO_STRAND, cmp_strand

# constants for fake 'start' and 'end' nodes
DUMMY_START_NODE = Exon(-1,-1)
DUMMY_END_NODE = Exon(-2,-2)

# constants for use as graph attributes
NODE_TSS_ID = 'tss_id'
NODE_WEIGHT = 'weight'

EDGE_DENSITY = 'density'
EDGE_OUT_FRAC = 'outfrac'
EDGE_IN_FRAC = 'infrac'

PLEN = 'plen'
PWEIGHT = 'pweight'
PDENSITY = 'pdensity'
PSRC = 'psrc'


def calculate_strand_fraction(G):
    '''
    compute the fraction of (+) versus (-) transcription in the
    graph by using the summed coverage densities in all the nodes 
    '''
    cluster_tree = ClusterTree(0,1)
    nodes = G.nodes()
    for i,n in enumerate(nodes):    
        cluster_tree.insert(n.start, n.end, i)
    for start, end, indexes in cluster_tree.getregions():
        strand_scores = [0, 0, 0]
        cluster_nodes = [nodes[i] for i in indexes]
        for n in cluster_nodes:
            for ndata in G.node[n]['data']:
                strand_scores[ndata.strand] += ndata.score / float((n.end - n.start))        
        total_score = strand_scores[POS_STRAND] + strand_scores[NEG_STRAND]        
        if total_score == 0:
            # if there is no "stranded" coverage at this node, then assign coverage
            # the positive strand by convention
            pos_frac = 1.0
        else:
            # proportionally assign unstranded coverage based on amount of
            # plus and minus strand coverage
            pos_frac = strand_scores[POS_STRAND] / float(total_score)
        #print 'CLUSTER_NODES', cluster_nodes, "FRAC", pos_frac
        yield cluster_nodes, pos_frac
    del cluster_tree

def sum_transcript_data_scores(tdata_list):
    '''
    sum exon data scores by strand
    '''
    # sum coverage by strand
    strand_weights = [0, 0, 0]
    for tdata in tdata_list:
        strand_weights[tdata.strand] += tdata.score
    return strand_weights

def build_strand_specific_graphs(G):
    '''
    build separate subgraphs of G - a forward strand graph and a 
    reverse strand graph (if necessary)
    
    pos_frac: fraction of all 'stranded' coverage that is on the 
    positive strand.  this is used to partition unstranded coverage
    proportionally
 
    nodes in the graph will have a 'weight' attribute equal to the
    total strand-specific coverage for that node

    edges (u,v) in the graph will have a 'weight' attribute corresponding
    to the fraction of total coverage flowing out of node 'u' that goes 
    to 'v'
    '''
    # make graphs
    GG = (nx.DiGraph(), nx.DiGraph())
    # cluster nodes and compute fraction of fwd/rev strand coverage    
    for cluster_nodes, pos_frac in calculate_strand_fraction(G):
        # add nodes first
        strand_fracs = (pos_frac, 1.0 - pos_frac)
        for n in cluster_nodes:
            d = G.node[n]
            # store strand fractions of this node
            d['strand_fracs'] = strand_fracs
            # partition the unstranded coverage to fwd/rev strands
            # according to the fraction of stranded coverage observed
            strand_weights = sum_transcript_data_scores(d['data'])
            for strand in (POS_STRAND, NEG_STRAND):
                w = strand_weights[strand] + strand_fracs[strand] * strand_weights[NO_STRAND] 
                if w > 0:
                    GG[strand].add_node(n, attr_dict={NODE_WEIGHT: w})        
            #print 'NODE', n, 'WEIGHTS', strand_weights, 'W', w                            
    # add edges
    for u,v,d in G.edges_iter(data=True):            
        # determine strand based on coordinates of nodes
        strand = int(u.start >= v.end)        
        # check against strand attribute
        # TODO: can remove this
        strands = set(x.strand for x in d['data'])
        assert all(cmp_strand(x, strand) for x in strands)
        # get strand fraction by averaging the strand fractions
        # of the source and destination nodes
        src_fracs = G.node[u]['strand_fracs']
        dst_fracs = G.node[v]['strand_fracs']
        strand_frac = (src_fracs[strand] + dst_fracs[strand]) / 2.0
        # partition the unstranded coverage to fwd/rev strands
        # according to the fraction of stranded coverage observed
        strand_densities = sum_transcript_data_scores(d['data'])
        density = strand_densities[strand] + (strand_frac * strand_densities[NO_STRAND]) 
        if density > 0:
            GG[strand].add_edge(u, v, attr_dict={EDGE_DENSITY: density})
        #print 'EDGE', u, v, 'STRAND', strand, 'WEIGHTS', strand_weights, 'W', w                            
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

def calculate_edge_attrs(G):
    for u,nbrdict in G.adjacency_iter():
        # find total coverage flowing out of this node
        out_density_total = sum(eattrs[EDGE_DENSITY] for eattrs in nbrdict.itervalues())        
        for v, eattrs in nbrdict.iteritems():
            # find fraction flowing out of each node
            out_frac = eattrs[EDGE_DENSITY]/float(out_density_total)
            eattrs[EDGE_OUT_FRAC] = out_frac            
    # now fraction of weight flowing in to each node
    # reverse edge directions and use adjacency iter again
    G.reverse(copy=False)
    for u,nbrdict in G.adjacency_iter():
        # find total coverage density flowing into this node
        in_density_total = sum(eattrs[EDGE_DENSITY] for eattrs in nbrdict.itervalues())
        for v, eattrs in nbrdict.iteritems():
            # find fraction flowing into node            
            in_frac = eattrs[EDGE_DENSITY]/float(in_density_total)
            eattrs[EDGE_IN_FRAC] = in_frac
    # reverse the edges back to normal
    G.reverse(copy=False)

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
    out_density_total = sum(G.edge[u][x][EDGE_DENSITY] for x in u_succs)
    if out_density_total == 0:
        out_density_total = 1.0
    # compute edge attributes leaving node 'u'
    for succ in u_succs:
        eattrs = G.edge[u][succ]
        out_frac = eattrs[EDGE_DENSITY]/float(out_density_total)
        eattrs[EDGE_OUT_FRAC] = out_frac
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
        # subtract density from edge
        u = path[i-1]
        eattrs = G.edge[u][v]
        eattrs[EDGE_DENSITY] = max(0, eattrs[EDGE_DENSITY] - path_density)
        # recompute edge attrs
        recalculate_edge_attrs(G, u, v)

def add_dummy_start_end_nodes(G, start_nodes, end_nodes):
    # add 'dummy' tss nodes if necessary
    G.add_node(DUMMY_START_NODE, weight=0)
    for start_node in start_nodes:        
        logging.debug('adding dummy %s -> %s' % (DUMMY_START_NODE, start_node))
        attr_dict = {EDGE_DENSITY: 0.0,
                     EDGE_OUT_FRAC: 0.0,
                     EDGE_IN_FRAC: 1.0}                     
        G.add_edge(DUMMY_START_NODE, start_node, attr_dict=attr_dict)
    # add a single 'dummy' end node
    G.add_node(DUMMY_END_NODE, weight=0)
    for end_node in end_nodes:
        logging.debug('adding dummy %s -> %s' % (end_node, DUMMY_END_NODE))
        attr_dict = {EDGE_DENSITY: 0.0,
                     EDGE_OUT_FRAC: 1.0,
                     EDGE_IN_FRAC: 0.0}
        G.add_edge(end_node, DUMMY_END_NODE, attr_dict=attr_dict)

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

        # DEBUG OUTPUT
        for u,nbrdict in G.adjacency_iter():
            logging.debug("NODE(%s) -> %f" % (u, G.node[u][NODE_WEIGHT]))
            for v, eattrs in nbrdict.iteritems():
                logging.debug("EDGE(%s,%s) -> %s" % (u, v, eattrs))
        
        # subtract the coverage of this path from the graph and recompute 
        # the edge weights
        subtract_path(G, path[1:-1], path_density)

        # DEBUG OUTPUT
        for u,nbrdict in G.adjacency_iter():
            logging.debug("SUBTRACTED NODE(%s) -> %f" % (u, G.node[u][NODE_WEIGHT]))
            for v, eattrs in nbrdict.iteritems():
                logging.debug("SUBTRACTED EDGE(%s,%s) -> %s" % (u, v, eattrs))
        
        
        
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
