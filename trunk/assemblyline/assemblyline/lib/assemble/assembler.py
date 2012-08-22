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

def add_dummy_nodes(G, dummy_nodes):
    # add to graph
    G.add_node(dummy_nodes[0], attr_dict={NODE_LENGTH:0, NODE_DENSITY:0})
    for i in xrange(1, len(dummy_nodes)):
        G.add_node(dummy_nodes[i], attr_dict={NODE_LENGTH:0, NODE_DENSITY:0})
        G.add_edge(dummy_nodes[i-1], dummy_nodes[i])

def add_dummy_start_end_nodes(G, k):
    """
    prepend/append 'k' dummy nodes to start/end of graph in order
    to facilitate the creation of a k-mer graph with a unique source
    and sink node where all paths are all least 'k' in length
    """
    # get all leaf nodes
    start_nodes = []
    end_nodes = []
    for n in G.nodes():
        if G.in_degree(n) == 0:
            start_nodes.append(n)
        if G.out_degree(n) == 0:
            end_nodes.append(n)   
    # create 'k' source dummy nodes
    source_node_id = -1
    source_dummy_nodes = [Exon(i,i) for i in xrange(source_node_id, source_node_id - k, -1)]
    # create 'k' sink dummy nodes
    sink_node_id = source_node_id - k
    sink_dummy_nodes = [Exon(i,i) for i in xrange(sink_node_id, sink_node_id - k, -1)]
    sink_dummy_nodes.reverse()
    # add dummy nodes to graph
    add_dummy_nodes(G, source_dummy_nodes)
    add_dummy_nodes(G, sink_dummy_nodes)
    # connect all leaf nodes to dummy nodes
    for n in start_nodes:
        G.add_edge(source_dummy_nodes[-1], n)
    for n in end_nodes:
        G.add_edge(n, sink_dummy_nodes[0])
    return tuple(source_dummy_nodes), tuple(sink_dummy_nodes)

def extend_partial_paths(G, partial_paths):
    """
    extend paths in both directions to include dummy start/end nodes
    """  
    new_partial_paths = []
    for path, density in partial_paths:
        new_path = tuple(path)
        # extend backward
        preds = G.predecessors(new_path[0])
        while (len(preds) == 1) and (preds[0].start < 0):
            new_path = (preds[0],) + new_path
            preds = G.predecessors(new_path[0])
        # extend forward
        succs = G.successors(new_path[-1])
        while (len(succs) == 1) and (succs[0].start < 0):
            new_path = new_path + (succs[0],)
            succs = G.successors(new_path[-1])
        new_partial_paths.append((new_path, density))
    return new_partial_paths

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

def hash_kmers(K, k, ksmall):
    kmerhash = collections.defaultdict(lambda: set())
    for kmer in K.nodes_iter():
        for i in xrange(k - (ksmall-1)):
            kmerhash[kmer[i:i+ksmall]].add(kmer)
    return kmerhash

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
        kmer_hash = hash_kmers(K, k, ksmall)
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

def assemble_transcript_graph(G, strand, partial_paths, kmax,
                              fraction_major_path, max_paths):
    """
    enumerates individual transcript isoforms from transcript graph using
    a greedy algorithm
    
    strand: strand of graph G
    fraction_major_path: only return isoforms with density greater than 
    some fraction of the highest density path
    max_paths: do not enumerate more than max_paths isoforms     
    """
    # don't allow fraction_major_path to equal zero because that 
    # will force each run to iterate to 'max_paths'.  also don't
    # allow to equal one because potentially no paths will be 
    # returned
    if fraction_major_path < 0:
        fraction_major_path = 0.0
    if fraction_major_path >= 1.0:
        fraction_major_path = 1.0
    # determine parameter 'k' for de bruijn graph assembly
    # kmax should not be set to less than '2'
    kmin = 2
    kmax = max(kmin, kmax)    
    # choose k >= kmin and k <= length of longest partial path
    longest_path = max(len(x[0]) for x in partial_paths)
    k = max(kmin, min(kmax, longest_path))
    logging.debug("\tConstructing k-mer graph with k=%d" % (k))
    # append/prepend dummy nodes to start/end nodes
    source_kmer, sink_kmer = add_dummy_start_end_nodes(G, k)
    # extend paths in both directions to include the dummy start/end nodes
    partial_paths = extend_partial_paths(G, partial_paths)   
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
