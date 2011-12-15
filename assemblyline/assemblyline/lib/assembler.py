'''
Created on Feb 20, 2011

@author: mkiyer
'''
import networkx as nx
import numpy as np

from transcript import Exon, NEG_STRAND, NO_STRAND
from collapse_chains import collapse_strand_specific_graph
from assembler_base import NODE_DENSITY, NODE_LENGTH, EDGE_OUT_FRAC, \
    EDGE_IN_FRAC, SOURCE_NODE, SINK_NODE, MIN_NODE_DENSITY, \
    GLOBAL_GENE_ID, GLOBAL_TSS_ID, GLOBAL_TRANSCRIPT_ID
from path_finder import find_suboptimal_paths, calculate_edge_fractions

class PathInfo(object):
    """object to store path finder results"""
    __slots__ = ("density", "tss_id", "tx_id", "path")
    def __init__(self, density, tss_id, path):
        self.density = density
        self.tss_id = tss_id
        self.path = path
        self.tx_id = -1

# private constant for use as graph attributes
NODE_TSS_ID = 'tss_id'

def find_start_and_end_nodes(G, strand):
    """
    finds start nodes and end nodes in graph and annotates start nodes with 
    'tss_id' attributes
    """
    global GLOBAL_TSS_ID
    
    tss_pos_id_map = {}
    start_nodes = []
    for n,d in G.in_degree_iter():
        # found TSS if no incoming edges
        if d == 0:
            tss_pos = n.end if strand == NEG_STRAND else n.start
            # map TSS positions to IDs
            if tss_pos not in tss_pos_id_map:
                tss_pos_id_map[tss_pos] = GLOBAL_TSS_ID
                GLOBAL_TSS_ID += 1
            # add to list of start nodes
            start_nodes.append(n)
            # add tss_id attribute to graph
            G.node[n][NODE_TSS_ID] = tss_pos_id_map[tss_pos]
    # find end nodes by checking for nodes with no outgoing edges
    end_nodes = [n for (n,d) in G.out_degree_iter()
                 if (d == 0)]
    return start_nodes, end_nodes

def add_dummy_start_end_nodes(G, start_nodes, end_nodes):
    # add 'dummy' tss nodes if necessary
    G.add_node(SOURCE_NODE, attr_dict={NODE_DENSITY: MIN_NODE_DENSITY, 
                                       NODE_LENGTH: 0})
    #source_attr_dict = {EDGE_IN_FRAC: 1.0,
    #                    EDGE_OUT_FRAC: 1.0/len(start_nodes)}
    for start_node in start_nodes:        
        G.add_edge(SOURCE_NODE, start_node)
        #logging.debug('adding dummy %s -> %s' % (SOURCE_NODE, start_node))
        #G.add_edge(SOURCE_NODE, start_node, attr_dict=source_attr_dict)
    # add a single 'dummy' end node
    G.add_node(SINK_NODE, attr_dict={NODE_DENSITY: MIN_NODE_DENSITY,
                                     NODE_LENGTH: 0})
    #sink_attr_dict = {EDGE_IN_FRAC: 1.0/len(start_nodes),
    #                  EDGE_OUT_FRAC: 1.0}
    for end_node in end_nodes:
        G.add_edge(end_node, SINK_NODE)
        #logging.debug('adding dummy %s -> %s' % (end_node, SINK_NODE))
        #G.add_edge(end_node, SINK_NODE, attr_dict=sink_attr_dict)

def expand_path_chains(G, strand, path):
    # reverse negative stranded data so that all paths go from 
    # small -> large genomic coords
    if strand == NEG_STRAND:
        path.reverse()
    # get chains (children nodes) along path
    newpath = []
    for n in path:
        newpath.extend(G.node[n]['chain'])
    path = newpath
    # collapse contiguous nodes along path
    newpath = []
    chain = [path[0]]
    for v in path[1:]:
        if chain[-1].end != v.start:
            # update path with merge chain node
            newpath.append(Exon(chain[0].start, 
                                chain[-1].end))
            # reset chain
            chain = []
        chain.append(v)
    # add last chain
    newpath.append(Exon(chain[0].start, 
                        chain[-1].end))
    return newpath

def assemble_subgraph(G, strand, fraction_major_path, max_paths):
    global GLOBAL_TRANSCRIPT_ID
    # find start and end nodes in graph and add the tss_id attribute
    # to all start nodes
    start_nodes, end_nodes = find_start_and_end_nodes(G, strand)
    # at an artificial node at the start and end so that all start
    # nodes are searched together for the best path
    add_dummy_start_end_nodes(G, start_nodes, end_nodes)
    # compute initial edge density and coverage flow in/out of nodes
    # calculate_edge_fractions(G, NODE_DENSITY)
    # TODO: redistribute node density in order to minimize error in graph
    # smoothen_graph_leastsq(G)
    # find up to 'max_paths' paths through graph
    path_info_list = []
    for path, density in find_suboptimal_paths(G, SOURCE_NODE, SINK_NODE, 
                                               fraction_major_path,
                                               max_paths):
        # remove dummy nodes from path
        path = list(path[1:-1])
        # get tss_id from path
        tss_id = G.node[path[0]][NODE_TSS_ID]
        # cleanup and expand chained nodes within path
        path = expand_path_chains(G, strand, path)
        # set transcript id
        path_info = PathInfo(density, tss_id, path)
        path_info.tx_id = GLOBAL_TRANSCRIPT_ID
        GLOBAL_TRANSCRIPT_ID += 1
        # add to path list
        path_info_list.append(path_info)
    # remove dummy nodes from graph
    G.remove_node(SOURCE_NODE)
    G.remove_node(SINK_NODE)
    return path_info_list

def assemble_transcript_graph(H, strand, fraction_major_path, max_paths):
    """
    returns (gene_id, list of PathInfo objects) tuple
    """
    global GLOBAL_GENE_ID
    # don't allow fraction_major_path to equal zero because that 
    # will force each run to iterate to 'max_paths'.  also don't
    # allow to equal one because potentially no paths will be 
    # returned
    if fraction_major_path <= 0:
        fraction_major_path = 1e-8
    if fraction_major_path >= 1.0:
        fraction_major_path = 1.0 - 1e-8
    # collapse chains of nodes in graph, which reformats the graph 
    # attributes, and then get connected components of graph which
    # represent independent genes
    G = collapse_strand_specific_graph(H)
    Gsubs = nx.weakly_connected_component_subgraphs(G)
    # get connected components - unconnected components are considered
    # different genes
    for Gsub in Gsubs:
        # get assembled paths
        path_info_list = assemble_subgraph(Gsub, strand,
                                           fraction_major_path=fraction_major_path, 
                                           max_paths=max_paths)
        yield GLOBAL_GENE_ID, path_info_list
        GLOBAL_GENE_ID += 1
