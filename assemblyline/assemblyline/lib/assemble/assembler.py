'''
Created on Feb 20, 2011

@author: mkiyer
'''
import logging
import collections
import networkx as nx

from assemblyline.lib.transcript import Exon, NEG_STRAND
from base import NODE_DENSITY, NODE_LENGTH, SMOOTH_FWD, SMOOTH_REV, SMOOTH_TMP, CHAIN_NODES
from path_finder import find_suboptimal_paths
from smooth import smooth_graph

class PathInfo(object):
    """object to store path finder results"""
    __slots__ = ("density", "path", "gene_id", "tss_id")
    def __init__(self, density, path):
        self.density = density
        self.path = path
        self.gene_id = -1
        self.tss_id = -1

def get_start_end_nodes(G):
    # get all leaf nodes
    start_nodes = set()
    end_nodes = set()
    for n in G.nodes_iter():
        if G.in_degree(n) == 0:
            start_nodes.add(n)
        if G.out_degree(n) == 0:
            end_nodes.add(n)
    return start_nodes, end_nodes

def get_dummy_nodes(k):
    # create 'k' source dummy nodes
    source_node_id = -1
    source_dummy_nodes = [Exon(i,i) for i in xrange(source_node_id, source_node_id - k, -1)]
    # create 'k' sink dummy nodes
    sink_node_id = source_node_id - k
    sink_dummy_nodes = [Exon(i,i) for i in xrange(sink_node_id, sink_node_id - k, -1)]
    sink_dummy_nodes.reverse()
    return tuple(source_dummy_nodes), tuple(sink_dummy_nodes)

def get_kmers(path, k):
    kmers = set()
    kmers.add(path[0:k])
    for i in xrange(0, len(path) - (k-1)):
        kmers.add(path[i:i+k])
    return kmers

def hash_kmers(kmers, k, ksmall):
    kmer_hash = collections.defaultdict(lambda: set())
    for kmer in kmers:
        for i in xrange(k - (ksmall-1)):
            kmer_hash[kmer[i:i+ksmall]].add(kmer)
    return kmer_hash

def optimize_k(G, partial_paths, kmax, sensitivity_threshold):
    """
    determine optimal choice for parameter 'k' for assembly
    maximizes k while ensuring sensitivity constraint is met
    """
    total_density = sum(density for path,density in partial_paths)
    total_paths = len(partial_paths)
    remaining_density = float(total_density)
    remaining_paths = list(partial_paths)
    k = 2
    # get all leaf nodes
    start_nodes, end_nodes = get_start_end_nodes(G)
    while k < kmax:
        # create 'k' dummy nodes
        source_dummy_nodes, sink_dummy_nodes = get_dummy_nodes(k)
        kmers = set()
        short_paths = []
        new_remaining_paths = []        
        for path, density in remaining_paths:
            # TODO: remove assert
            assert len(path) > 0
            # extend paths that start or end on source or sink nodes
            new_path = tuple(path)
            if path[0] in start_nodes:
                new_path = source_dummy_nodes + new_path
            if path[-1] in end_nodes:
                new_path = new_path + sink_dummy_nodes
            if len(new_path) < k:
                short_paths.append((new_path, density))
            else:
                kmers.update(get_kmers(new_path,k))
                new_remaining_paths.append((path,density))
        # create a (k-1) to k-mer hash
        kmer_hash = hash_kmers(kmers, k, k-1)
        for path,density in short_paths:
            if path not in kmer_hash:
                remaining_density -= density
                continue
        sensitivity = (remaining_density / total_density)
        logging.debug("\t\tk=%d sensitivity=%f (%d/%d paths) (%f/%f density)" %
                      (k, sensitivity, len(new_remaining_paths), 
                       total_paths, remaining_density, total_density))
        if sensitivity < sensitivity_threshold:
            break
        k += 1
        remaining_paths = new_remaining_paths
    return k

def add_dummy_nodes(G, dummy_nodes):
    # add to graph
    G.add_node(dummy_nodes[0], attr_dict={NODE_LENGTH:0, NODE_DENSITY:0})
    for i in xrange(1, len(dummy_nodes)):
        G.add_node(dummy_nodes[i], attr_dict={NODE_LENGTH:0, NODE_DENSITY:0})
        G.add_edge(dummy_nodes[i-1], dummy_nodes[i])

def attach_dummy_leaf_nodes(G, k, partial_paths):
    """
    prepend/append 'k' dummy nodes to start/end of graph in order
    to facilitate the creation of a k-mer graph with a unique source
    and sink node where all paths are all least 'k' in length

    also extend paths in both directions to include dummy start/end nodes
    """
    # get all leaf nodes
    start_nodes, end_nodes = get_start_end_nodes(G)
    # create 'k' dummy nodes
    source_dummy_nodes, sink_dummy_nodes = get_dummy_nodes(k)    
    # extend paths that start or end on source or sink nodes
    new_partial_paths = []
    for path, density in partial_paths:
        new_path = tuple(path)
        if path[0] in start_nodes:
            new_path = source_dummy_nodes + new_path
        if path[-1] in end_nodes:
            new_path = new_path + sink_dummy_nodes
        new_partial_paths.append((new_path, density))
    # add dummy nodes to graph
    add_dummy_nodes(G, source_dummy_nodes)
    add_dummy_nodes(G, sink_dummy_nodes)
    # connect all leaf nodes to dummy nodes
    for n in start_nodes:
        G.add_edge(source_dummy_nodes[-1], n)
    for n in end_nodes:
        G.add_edge(n, sink_dummy_nodes[0])
    return tuple(source_dummy_nodes), tuple(sink_dummy_nodes), new_partial_paths

def init_kmer_attrs():
    return {NODE_DENSITY: 0.0, SMOOTH_FWD: 0.0, 
            SMOOTH_REV: 0.0, SMOOTH_TMP: 0.0}
    
def add_path(K, k, path, density):
    # break longer paths into kmers of length 'k' and add them
    # to the graph
    from_kmer = path[0:k]
    if from_kmer not in K:
        K.add_node(from_kmer, attr_dict=init_kmer_attrs())
    kmerattrs = K.node[from_kmer]
    kmerattrs[NODE_DENSITY] += density
    # the first kmer should be "smoothed" in reverse direction
    kmerattrs[SMOOTH_REV] += density
    for i in xrange(1,len(path) - (k-1)):
        to_kmer = path[i:i+k]
        if to_kmer not in K:
            K.add_node(to_kmer, attr_dict=init_kmer_attrs())
        kmerattrs = K.node[to_kmer]
        kmerattrs[NODE_DENSITY] += density
        # connect kmers
        K.add_edge(from_kmer, to_kmer)
        # update from_kmer to continue loop
        from_kmer = to_kmer
    # the last kmer should be "smoothed" in forward direction
    kmerattrs[SMOOTH_FWD] += density

def extrapolate_short_path(kmer_hash, kmer_density_dict, path, density):
    """
    add the list of 'partial_paths' paths to graph 'K'

    paths shorter than 'k' are extrapolated into k-mers before being added 
    to the graph
    """
    if path not in kmer_hash:
        return []
    total_density = 0
    matching_kmers = []
    for kmer in kmer_hash[path]:
        # get approximate density at kmer
        kmer_density = kmer_density_dict[kmer]        
        # compute total density at matching kmers
        total_density += kmer_density
        matching_kmers.append((kmer, kmer_density))
    # now calculate fractional densities for matching kmers
    new_partial_paths = [] 
    for kmer,kmer_density in matching_kmers:
        new_density = density * (kmer_density / total_density)
        new_partial_paths.append((kmer, new_density))
    return new_partial_paths

def create_kmer_graph(G, partial_paths, k):
    """
    add partial paths from graph 'G' to k-mer graph 'K'

    partial_paths is a list of (path,density) tuples
    """
    # create a graph of k-mers
    K = nx.DiGraph()
    short_partial_path_dict = collections.defaultdict(lambda: [])
    for path, density in partial_paths:
        if len(path) < k:
            short_partial_path_dict[len(path)].append((path, density))
        else:
            add_path(K, k, path, density)
    # get kmer densities from linear graph
    kmer_density_dict = {}
    for kmer in K.nodes_iter():
        node_densities = [G.node[n][NODE_DENSITY] for n in kmer if n.start >= 0]
        if len(node_densities) > 0:
            kmer_density_dict[kmer] = min(node_densities)
    # extrapolate short kmers and add to graph
    for ksmall, short_partial_paths in short_partial_path_dict.iteritems():
        kmer_hash = hash_kmers(K.nodes_iter(), k, ksmall)
        for path, density in short_partial_paths:
            new_paths = extrapolate_short_path(kmer_hash, kmer_density_dict, path, density)
            for p,d in new_paths:
                add_path(K, k, p, d)
    return K

def connect_source_and_sink_kmers(K, source_kmer, sink_kmer):
    # connect all kmer nodes with degree zero to the source/sink node
    # to account for fragmentation in the kmer graph when k > 2
    k = len(source_kmer)
    kmers = K.nodes()
    for kmer in kmers:
        if K.in_degree(kmer) == 0:
            if kmer == source_kmer:
                continue
            add_path(K, k, source_kmer + kmer, 0.0)
        if K.out_degree(kmer) == 0:
            if kmer == sink_kmer:
                continue
            add_path(K, k, kmer + sink_kmer, 0.0)

def expand_path_chains(G, strand, path):
    # reverse negative stranded data so that all paths go from 
    # small -> large genomic coords
    if strand == NEG_STRAND:
        path.reverse()
    # get chains (children nodes) along path
    newpath = []
    for n in path:
        nd = G.node[n]
        if CHAIN_NODES in nd:
            newpath.extend(nd[CHAIN_NODES])
        else:
            newpath.append(n)
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
    newpath.append(Exon(chain[0].start, chain[-1].end))
    return newpath

def assemble_transcript_graph(G, strand, partial_paths, 
                              sensitivity_threshold, user_kmax,
                              fraction_major_path, 
                              max_paths):
    """
    enumerates individual transcript isoforms from transcript graph using
    a greedy algorithm
    
    strand: strand of graph G
    fraction_major_path: only return isoforms with density greater than 
    some fraction of the highest density path
    max_paths: do not enumerate more than max_paths isoforms     
    """
    # constrain sensitivity parameter
    sensitivity_threshold = min(max(0.0, sensitivity_threshold), 1.0)
    # constrain fraction_major_path parameter
    fraction_major_path = min(max(0.0, fraction_major_path), 1.0)
    # determine parameter 'k' for de bruijn graph assembly
    kmax = max(len(x[0]) for x in partial_paths)
    if user_kmax is not None:
        kmax = min(kmax, user_kmax)
    logging.debug("\tOptimizing k parameter (kmax=%d)" % (kmax))
    k = optimize_k(G, partial_paths, kmax, sensitivity_threshold)
    logging.debug("\tConstructing k-mer graph with k=%d" % (k))
    # append/prepend dummy nodes to start/end nodes
    source_kmer, sink_kmer, partial_paths = \
        attach_dummy_leaf_nodes(G, k, partial_paths)
    # create k-mer graph and add partial paths
    K = create_kmer_graph(G, partial_paths, k)
    # connect all kmer nodes with degree zero to the source/sink node
    # to account for fragmentation in the kmer graph when k > 2
    connect_source_and_sink_kmers(K, source_kmer, sink_kmer)
    # smooth kmer graph
    smooth_graph(K)
    # find up to 'max_paths' paths through graph
    logging.debug("\tFinding suboptimal paths in k-mer graph with %d nodes" % (len(K)))
    path_info_list = []   
    for kmer_path, density in find_suboptimal_paths(K, source_kmer, 
                                                    sink_kmer,
                                                    fraction_major_path, 
                                                    max_paths):
        # get nodes from kmer path
        path = list(kmer_path[0])
        path.extend(kmer[-1] for kmer in kmer_path[1:])
        # remove dummy nodes
        path = [n for n in path if n.start >= 0]
        # cleanup and expand nodes within path
        path = expand_path_chains(G, strand, path)
        # add to path list
        path_info_list.append(PathInfo(density, path))
        logging.debug("\t\tdensity=%f num_exons=%d" % (density, len(path)))
    return path_info_list
