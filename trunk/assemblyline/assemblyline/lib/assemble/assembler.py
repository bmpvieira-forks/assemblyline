'''
Created on Feb 20, 2011

@author: mkiyer
'''
import logging
import collections
import networkx as nx

from assemblyline.lib.transcript import Exon, NEG_STRAND
from base import NODE_DENSITY, NODE_LENGTH 
from path_finder import find_suboptimal_paths

# private constants
SMOOTH_FWD = 'smfwd'
SMOOTH_REV = 'smrev'
SMOOTH_TMP = 'smtmp'

class PathInfo(object):
    """object to store path finder results"""
    __slots__ = ("density", "path", "gene_id", "tss_id")
    def __init__(self, density, path):
        self.density = density
        self.path = path
        self.gene_id = -1
        self.tss_id = -1

def _init_kmer_attrs():
    return {NODE_DENSITY: 0.0, SMOOTH_FWD: 0.0, 
            SMOOTH_REV: 0.0, SMOOTH_TMP: 0.0}

def _extend_partial_paths(G, partial_paths):
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

def _add_dummy_nodes(G, dummy_nodes):
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
    _add_dummy_nodes(G, source_dummy_nodes)
    _add_dummy_nodes(G, sink_dummy_nodes)
    # connect all leaf nodes to dummy nodes
    for n in start_nodes:
        G.add_edge(source_dummy_nodes[-1], n)
    for n in end_nodes:
        G.add_edge(n, sink_dummy_nodes[0])
    return tuple(source_dummy_nodes), tuple(sink_dummy_nodes)

def get_forward_kmers(G, seed_path, k):
    paths = []
    stack = [seed_path]
    while len(stack) > 0:
        path = stack.pop()
        if len(path) == k:
            paths.append(path)
        else:
            for succ in G.successors_iter(path[-1]):
                stack.append(path + (succ,))
    return paths

#def old_create_kmer_graph(G, k):
#    """
#    returns a DiGraph of all possible k-mers    
#    """
#    # create a graph of k-mers
#    K = nx.DiGraph()    
#    # TODO: use integer nodes instead of Exon objects to improve 
#    # performance
#    kmer_hash = collections.defaultdict(lambda: (set(),set()))   
#    for n in G.nodes_iter():
#        # get all kmers that start with node
#        kmers = get_forward_kmers(G, (n,), k)
#        # add kmers to (k-1) hash
#        for kmer in kmers:
#            # add node
#            K.add_node(kmer, attr_dict={NODE_DENSITY: 0.0,
#                                        SMOOTH_FWD: 0.0,
#                                        SMOOTH_REV: 0.0,
#                                        SMOOTH_TMP: 0.0})
#            # update hash
#            kmer_hash[kmer[1:]][REV].add(kmer[0])
#            kmer_hash[kmer[:-1]][FWD].add(kmer[-1])
#    # now connect kmers using kmer hash
#    for kminus1mer, fwd_rev_node_sets in kmer_hash.iteritems():
#        for fwd_node in fwd_rev_node_sets[FWD]:
#            for rev_node in fwd_rev_node_sets[REV]:
#                # add edge
#                from_kmer = (rev_node,) + kminus1mer 
#                to_kmer = kminus1mer + (fwd_node,)
#                K.add_edge(from_kmer, to_kmer)
#    return K

def _bin_partial_paths(G, partial_paths, k):
    """
    bin partial paths by their length (up to length 'k') and sort 
    paths in each bin by their density levels
    """   
    # bin paths by length
    path_bins = collections.defaultdict(lambda: [])
    for path, density in partial_paths:
        kbin = k if len(path) > k else len(path)
        path_bins[kbin].append((path,density))
    # sort paths within each bin by density
    for kbin in path_bins:
        path_bins[kbin].sort(key=lambda x: x[1])
    return path_bins

def _add_path(K, path, density, k):
    # break longer paths into kmers of length 'k' and add them
    # to the graph
    from_kmer = path[0:k]
    if from_kmer not in K:
        K.add_node(from_kmer, attr_dict=_init_kmer_attrs())
    kmerattrs = K.node[from_kmer]
    kmerattrs[NODE_DENSITY] += density
    # the first kmer should be "smoothed" in reverse direction
    kmerattrs[SMOOTH_REV] += density
    for i in xrange(1,len(path) - (k-1)):
        to_kmer = path[i:i+k]
        if to_kmer not in K:
            K.add_node(to_kmer, attr_dict=_init_kmer_attrs())
        kmerattrs = K.node[to_kmer]
        kmerattrs[NODE_DENSITY] += density
        # connect kmers
        K.add_edge(from_kmer, to_kmer)
        # update from_kmer to continue loop
        from_kmer = to_kmer
    # the last kmer should be "smoothed" in forward direction
    kmerattrs[SMOOTH_FWD] += density

def _add_long_partial_paths(G, K, partial_paths, k):
    """
    add the list of 'partial_paths' paths to graph 'K'

    paths longer than 'k' are broken into k-mers before being 
    added to the graph. 
    """
    while len(partial_paths) > 0:
        # get next highest density path from list
        path, density = partial_paths[-1]
        _add_path(K, path, density, k)
        # pop path from list
        partial_paths.pop()

def _extend_path_reverse(G, seed, length):
    paths = []
    stack = [(seed,)]
    while len(stack) > 0:
        path = stack.pop()
        if len(path) == (length+1):
            paths.append(path[:-1])
        else:
            for pred in G.predecessors_iter(path[0]):
                stack.append((pred,) + path)
    return paths

def _extend_path_forward(G, seed, length):
    paths = []
    stack = [(seed,)]
    while len(stack) > 0:
        path = stack.pop()
        if len(path) == (length+1):
            paths.append(path[1:])
        else:
            for succ in G.successors_iter(path[-1]):
                stack.append(path + (succ,))
    return paths

def _get_partial_path_kmers(G, seed_path, k):
    assert len(seed_path) < k
    # since seed_path length is less than k, need to extend to a set
    # of possible k-mers
    L = (k - len(seed_path))
    Lrev = (L+1) >> 1
    Lfwd = L >> 1
    # extend path in both directions
    rev_paths = _extend_path_reverse(G, seed_path[0], Lrev)
    fwd_paths = _extend_path_forward(G, seed_path[-1], Lfwd)
    # enumerate kmers
    kmers = []
    for rp in rev_paths:
        for fp in fwd_paths:
            kmers.append(rp + seed_path + fp)
    return kmers

def _add_short_partial_paths(G, K, partial_paths, k):
    """
    add the list of 'partial_paths' paths to graph 'K'

    paths shorter than 'k' are extrapolated into k-mers before being added 
    to the graph
    """
    while len(partial_paths) > 0:
        # get next highest density path from list
        path, density = partial_paths[-1]
        # extend path in both directions to find all compatible 'kmers'
        kmers = _get_partial_path_kmers(G, path, k)
        # compute total density at matching kmers
        total_density = 0
        new_kmers = []
        for kmer in kmers:
            # restrict to kmers that have already been added to graph K
            if kmer not in K:
                continue
            total_density += K.node[kmer][NODE_DENSITY]
            new_kmers.append(kmer)
        # apply density proportionately to all matching kmers 
        for kmer in new_kmers:
            kmerattrs = K.node[kmer]
            frac = kmerattrs[NODE_DENSITY] / float(total_density)
            adj_density = (frac * density)
            kmerattrs[NODE_DENSITY] += adj_density
            kmerattrs[SMOOTH_FWD] += adj_density
            kmerattrs[SMOOTH_REV] += adj_density
        # pop path from list
        partial_paths.pop()

def create_kmer_graph(G, partial_paths, k):
    """
    add partial paths from graph 'G' to k-mer graph 'K'

    partial_paths is a list of (path,density) tuples
    """
    # create a graph of k-mers
    K = nx.DiGraph()
    # bin partial paths by their length (up to length 'k') and sort 
    # paths in each bin by their density levels
    partial_path_bins = _bin_partial_paths(G, partial_paths, k)
    # add partial paths of length >= k
    _add_long_partial_paths(G, K, partial_path_bins[k], k)
    del partial_path_bins[k]
    # iterate from large to small path-length bins
    for kbin in sorted(partial_path_bins, reverse=True):
        _add_short_partial_paths(G, K, partial_path_bins[kbin], k)
    return K

def _smooth_iteration(K, density_attr, smooth_attr):
    nodes = nx.topological_sort(K)
    for u in nodes:
        ud = K.node[u]
        smooth_density = ud[smooth_attr]
        succ = K.successors(u)
        if len(succ) == 0:
            continue
        total_nbr_density = sum(K.node[v][density_attr] for v in succ)        
        if total_nbr_density == 0:
            # if all successors have zero density apply smoothing evenly
            avg_density = smooth_density / len(succ)
            for v in succ:
                vd = K.node[v]
                vd[SMOOTH_TMP] += avg_density
                vd[smooth_attr] += avg_density
        else:
            # apply smoothing proportionately
            for v in succ:
                vd = K.node[v]
                frac = vd[density_attr]/float(total_nbr_density)
                adj_density = frac * smooth_density
                vd[SMOOTH_TMP] += adj_density
                vd[smooth_attr] += adj_density

def smooth_kmer_graph(K, density_attr=NODE_DENSITY):
    # smooth in forward direction
    _smooth_iteration(K, density_attr, SMOOTH_FWD)
    # smooth in reverse direction
    K.reverse(copy=False)
    _smooth_iteration(K, density_attr, SMOOTH_REV)
    K.reverse(copy=False)
    # apply densities to nodes
    for n,d in K.nodes_iter(data=True):
        d[density_attr] += d[SMOOTH_TMP]
   
def connect_source_and_sink_kmers(K, source_kmer, sink_kmer):
    # connect all kmer nodes with degree zero to the source/sink node
    # to account for fragmentation in the kmer graph when k > 2
    k = len(source_kmer)
    kmers = K.nodes()
    for kmer in kmers:
        if K.in_degree(kmer) == 0:
            if kmer == source_kmer:
                continue
            _add_path(K, source_kmer + kmer, 0.0, k)
        if K.out_degree(kmer) == 0:
            if kmer == sink_kmer:
                continue
            _add_path(K, kmer + sink_kmer, 0.0, k)

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
    newpath.append(Exon(chain[0].start, chain[-1].end))
    return newpath

def assemble_transcript_graph(G, strand, partial_paths, kmax, fraction_major_path, max_paths):
    """
    enumerates individual transcript isoforms from transcript graph using
    a greedy algorithm
    
    strand: strand of graph G
    partial_paths: list of (path,density) tuples comprising transcripts
    kmax: de-bruijn graph assembly parameter
    fraction_major_path: only return isoforms with density greater than 
    some fraction of the highest density path
    max_paths: do not enumerate more than max_paths isoforms     
    """
    global GLOBAL_TRANSCRIPT_ID
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
    partial_paths = _extend_partial_paths(G, partial_paths)   
    # create k-mer graph and add partial paths
    K = create_kmer_graph(G, partial_paths, k)
    # connect all kmer nodes with degree zero to the source/sink node
    # to account for fragmentation in the kmer graph when k > 2
    connect_source_and_sink_kmers(K, source_kmer, sink_kmer)
    # smooth kmer graph
    smooth_kmer_graph(K)
    # find up to 'max_paths' paths through graph
    logging.debug("\tFinding suboptimal paths in k-mer graph with %d nodes" % (len(K)))
    path_info_list = []
    for kmer_path, density in find_suboptimal_paths(K, fraction_major_path, 
                                                    max_paths):
        # get nodes from kmer path
        path = list(kmer_path[0])
        path.extend(kmer[-1] for kmer in kmer_path[1:])
        # remove dummy nodes
        path = [n for n in path if n.start >= 0]
        # cleanup and expand chained nodes within path
        path = expand_path_chains(G, strand, path)
        # add to path list
        path_info_list.append(PathInfo(density, path))
        logging.debug("\t\tdensity=%f length=%d nodes" % (density, len(path)))
    return path_info_list
