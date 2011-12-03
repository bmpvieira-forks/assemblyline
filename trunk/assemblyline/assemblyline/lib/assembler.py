'''
Created on Feb 20, 2011

@author: mkiyer
'''
import collections
import logging
import operator
import networkx as nx
import numpy as np

from bx.cluster import ClusterTree

from transcript import Exon, POS_STRAND, NEG_STRAND, NO_STRAND, cmp_strand
from smoothen import smoothen_graph_leastsq
from path_finder import find_suboptimal_paths

# assembler graph attributes
NODE_WEIGHT = 'weight'
EDGE_OUT_FRAC = 'outfrac'
EDGE_IN_FRAC = 'infrac'
EDGE_DENSITY = 'density'

# constants for fake 'start' and 'end' nodes
SOURCE_NODE = Exon(-1,-1)
SINK_NODE = Exon(-2,-2)

# private constants for use as graph attributes
NODE_TSS_ID = 'tss_id'
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
            # if there is no "stranded" coverage at this node, 
            # then assign coverage the positive strand by convention
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
            #print 'NODE', n, 'FRACS', strand_fracs, 'WEIGHTS', strand_weights, 'W', w                            
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

def add_dummy_start_end_nodes(G, start_nodes, end_nodes):
    # add 'dummy' tss nodes if necessary
    G.add_node(SOURCE_NODE, weight=0)
    for start_node in start_nodes:        
        logging.debug('adding dummy %s -> %s' % (SOURCE_NODE, start_node))
        attr_dict = {EDGE_DENSITY: 0.0,
                     EDGE_OUT_FRAC: 0.0,
                     EDGE_IN_FRAC: 1.0}                     
        G.add_edge(SOURCE_NODE, start_node, attr_dict=attr_dict)
    # add a single 'dummy' end node
    G.add_node(SINK_NODE, weight=0)
    for end_node in end_nodes:
        logging.debug('adding dummy %s -> %s' % (end_node, SINK_NODE))
        attr_dict = {EDGE_DENSITY: 0.0,
                     EDGE_OUT_FRAC: 1.0,
                     EDGE_IN_FRAC: 0.0}
        G.add_edge(end_node, SINK_NODE, attr_dict=attr_dict)

def smoothen_graph_leastsq(G):
    # cannot smooth a graph with just one node
    if len(G) == 1:
        return
    # label nodes
    cur_id = 0
    eqns = {}
    cols = {}
    y = np.zeros(len(G), dtype=np.float)
    # assign numeric ids to nodes
    for n,d in G.nodes_iter(data=True):
        d['id'] = cur_id
        y[cur_id] = d[NODE_WEIGHT]
        if G.in_degree(n) > 0:
            a = np.zeros(len(G), dtype=np.float)
            a[cur_id] = -1.0
            eqns[cur_id] = a
        cols[cur_id] = n
        cur_id += 1
    # setup equations for each node
    for u,v,d in G.edges_iter(data=True):
        u_id = G.node[u]['id']
        v_id = G.node[v]['id']
        out_frac = d[EDGE_OUT_FRAC]
        eqns[v_id][u_id] = out_frac
    # build matrix for least squares error
    # minimization
    row_ids = sorted(eqns)
    A = np.vstack(eqns[id] for id in row_ids)
    A = np.asmatrix(A)
    I = np.identity(len(G), dtype=np.float)
    # perform least squares minimization
    beta = (I - A.T*((A*A.T).I)*A).dot(y)
    # set new node weights
    for col_id,n in cols.iteritems():
        G.node[n][NODE_WEIGHT] = beta[0,col_id]
        del G.node[n]['id']
    #for n,d in G.nodes_iter(data=True):
    #    print 'NEW NODE', n, d

def assemble_subgraph(G, strand, fraction_major_path, 
                      max_paths):
    # find start and end nodes in graph
    start_nodes, end_nodes = find_start_and_end_nodes(G, strand)
    logging.debug("START NODES (TSSs): %s" % start_nodes)
    logging.debug("END NODES: %s" % end_nodes)
    # compute initial edge weights and coverage flow in/out of nodes
    calculate_edge_attrs(G)
    # redistribute node weights in order to minimize error in graph
    smoothen_graph_leastsq(G)
    # at an artificial node at the start and end so that all start
    # nodes are searched together for the best path
    add_dummy_start_end_nodes(G, start_nodes, end_nodes)
    # find up to 'max_paths' paths through graph
    paths = []
    for path, weight, length in \
        find_suboptimal_paths(G, SOURCE_NODE, SINK_NODE, max_paths):
        # remove dummy nodes from path
        path = path[1:-1]
        # compute density
        density = 1.0e3 * weight / float(length)
        # get tss_id from path
        tss_id = G.node[path[0]][NODE_TSS_ID]
        paths.append((density, tss_id, path))
    # sort by density
    paths = sorted(paths, key=operator.itemgetter(0), reverse=True)
    # return paths while density is greater than some fraction of the
    # highest density
    density_lower_bound = fraction_major_path * paths[0][0]
    for density, tss_id, path in paths:
        if density < density_lower_bound:
            break
        yield tss_id, density, path
    # remove dummy nodes from graph
    G.remove_node(SOURCE_NODE)
    G.remove_node(SINK_NODE)


def assemble_transcript_graph(G, fraction_major_path, max_paths):
    if fraction_major_path <= 0:
        fraction_major_path = 1e-8        
    # transform transcript graph into strand-specific graphs
    # for forward and reverse strands, and calculate node and
    # edge weights
    GG = build_strand_specific_graphs(G)
    # assemble strands one at a time
    current_gene_id = 0
    current_tss_id = 0 
    for strand, Gstrand in enumerate(GG):
        logging.debug("STRAND: %s" % strand) 
        # get connected components - unconnected components are considered
        # different genes
        for Gsub in nx.weakly_connected_component_subgraphs(Gstrand):
            tss_id_map = {}
            for sub_tss_id, score, path in \
                assemble_subgraph(Gsub, strand,
                                  fraction_major_path=fraction_major_path, 
                                  max_paths=max_paths):
                if sub_tss_id not in tss_id_map:
                    tss_id_map[sub_tss_id] = current_tss_id
                    current_tss_id += 1
                tss_id = tss_id_map[sub_tss_id]
                yield strand, current_gene_id, tss_id, score, path
            current_gene_id += 1

