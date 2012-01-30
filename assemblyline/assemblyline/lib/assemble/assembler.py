'''
Created on Feb 20, 2011

@author: mkiyer
'''
import networkx as nx
import collections
import numpy as np

from assemblyline.lib.transcript import Exon, NEG_STRAND, NO_STRAND
from base import NODE_DENSITY, NODE_LENGTH, NODE_TSS_ID, FWD, REV, EDGE_IN_FRAC, EDGE_OUT_FRAC, \
    GLOBAL_GENE_ID, GLOBAL_TSS_ID, GLOBAL_TRANSCRIPT_ID
from path_finder import find_suboptimal_paths

# private constants
DUMMY_NODE_ID = 1
SMOOTH_FWD = 'smfwd'
SMOOTH_REV = 'smrev'
SMOOTH_TMP = 'smtmp'

class PathInfo(object):
    """object to store path finder results"""
    __slots__ = ("density", "tss_id", "tx_id", "path")
    def __init__(self, density, tss_id, path):
        self.density = density
        self.tss_id = tss_id
        self.path = path
        self.tx_id = -1

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

def add_dummy_nodes(G, n, k, dir=FWD):
    global DUMMY_NODE_ID
    if k <= 1:
        return n
    # create 'k-1' dummy nodes
    dummy_nodes = [Exon(-id,-id) for id in xrange(DUMMY_NODE_ID, DUMMY_NODE_ID + (k-1))]
    DUMMY_NODE_ID += (k-1)
    # add to graph
    G.add_node(dummy_nodes[0], attr_dict={NODE_LENGTH:0, NODE_DENSITY:0})
    for i in xrange(1,k-1):
        G.add_node(dummy_nodes[i], attr_dict={NODE_LENGTH:0, NODE_DENSITY:0})
        G.add_edge(dummy_nodes[i-1], dummy_nodes[i])
    # link to true node
    if dir == FWD:
        G.add_edge(n, dummy_nodes[0])
        # return end node
        return dummy_nodes[-1]
    else:
        # return start node
        G.add_edge(dummy_nodes[-1], n)
        return dummy_nodes[0]

def add_dummy_start_end_nodes(G, k):
    """
    returns lists of start nodes and end nodes
    """
    start_nodes = set()
    end_nodes = set()
    nodes = G.nodes()
    for n in nodes:
        if G.in_degree(n) == 0:
            start_nodes.add(add_dummy_nodes(G, n, k, REV))
        if G.out_degree(n) == 0:
            end_nodes.add(add_dummy_nodes(G, n, k, FWD))
    return start_nodes, end_nodes

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

def create_kmer_graph(G, k):
    """
    returns a DiGraph of all possible k-mers    
    """
    # create a graph of k-mers
    K = nx.DiGraph()    
    # TODO: use integer nodes instead of Exon objects to improve 
    # performance
    kmer_hash = collections.defaultdict(lambda: (set(),set()))   
    for n in G.nodes_iter():
        # get all kmers that start with node
        kmers = get_forward_kmers(G, (n,), k)
        # add kmers to (k-1) hash
        for kmer in kmers:
            # add node
            K.add_node(kmer, attr_dict={NODE_DENSITY: 0.0,
                                        SMOOTH_FWD: 0.0,
                                        SMOOTH_REV: 0.0,
                                        SMOOTH_TMP: 0.0})
            # update hash (not necessary if k equals 1)
            if k > 1:
                #kmer[0] -> kmer[1:]
                kmer_hash[kmer[1:]][REV].add(kmer[0])
                #kmer[:-1] -> kmer[-1]
                kmer_hash[kmer[:-1]][FWD].add(kmer[-1])
    # now connect kmers using kmer hash
    for kminus1mer, fwd_rev_node_sets in kmer_hash.iteritems():
        for fwd_node in fwd_rev_node_sets[FWD]:
            for rev_node in fwd_rev_node_sets[REV]:
                # add edge
                from_kmer = (rev_node,) + kminus1mer 
                to_kmer = kminus1mer + (fwd_node,)
                K.add_edge(from_kmer, to_kmer)
    return K

def extend_path_reverse(G, seed, length):
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

def extend_path_forward(G, seed, length):
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

def get_partial_path_kmers(G, seed_path, k):
    assert len(seed_path) < k
    # since seed_path length is less than k, need to extend to a set
    # of possible k-mers
    L = (k - len(seed_path))
    Lrev = (L+1) >> 1
    Lfwd = L >> 1
    # extend path in both directions
    rev_paths = extend_path_reverse(G, seed_path[0], Lrev)
    fwd_paths = extend_path_forward(G, seed_path[-1], Lfwd)
    # enumerate kmers
    kmers = []
    for rp in rev_paths:
        for fp in fwd_paths:
            kmers.append(rp + seed_path + fp)
    return kmers

def add_partial_path_kmers(G, K, partial_paths, k):
    """
    add partial paths from graph 'G' to 'k'-mer graph 'K'

    partial_paths must be sorted by smallest -> largest path fragment
    """
    # start by adding partial paths with length >= k
    while len(partial_paths) > 0:
        # get next largest partial path from list
        path, density = partial_paths[-1]
        # stop iteration when path lengths drop below 'k'
        if len(path) < k:
            break
        # create kmers and add to graph
        kmer = path[0:k]
        assert kmer in K
        kmerattrs = K.node[kmer]
        kmerattrs[NODE_DENSITY] += density
        # the first kmer should be "smoothed" in reverse direction
        kmerattrs[SMOOTH_REV] += density
        for i in xrange(1,len(path) - (k-1)):
            kmer = path[i:i+k]
            assert kmer in K
            kmerattrs = K.node[kmer]
            kmerattrs[NODE_DENSITY] += density
        # the last kmers should be "smoothed"
        kmerattrs[SMOOTH_FWD] += density
        # pop path from list
        partial_paths.pop()
    # now add partial paths with length < k
    while len(partial_paths) > 0:
        # get next largest partial path from list
        path, density = partial_paths.pop()
        # extend short partial paths in both directions to find all 
        # compatible kmers
        kmers = get_partial_path_kmers(G, path, k)
        total_density = sum(K.node[kmer][NODE_DENSITY] for kmer in kmers)
        if total_density == 0:
            # all k-mers have zero density, so divide the density equally
            avg_density = density / float(len(kmers))
            for kmer in kmers:
                kmerattrs = K.node[kmer]
                kmerattrs[NODE_DENSITY] += avg_density
                kmerattrs[SMOOTH_FWD] += avg_density
                kmerattrs[SMOOTH_REV] += avg_density
        else:
            # apply density proportionately to all matching kmers 
            for kmer in kmers:
                kmerattrs = K.node[kmer]
                frac = kmerattrs[NODE_DENSITY] / float(total_density)
                adj_density = (frac * density)
                kmerattrs[NODE_DENSITY] += adj_density
                kmerattrs[SMOOTH_FWD] += adj_density
                kmerattrs[SMOOTH_REV] += adj_density

def smooth_iteration(K, density_attr, smooth_attr):
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
    smooth_iteration(K, density_attr, SMOOTH_FWD)
    # smooth in reverse direction
    K.reverse(copy=False)
    smooth_iteration(K, density_attr, SMOOTH_REV)
    K.reverse(copy=False)
    # apply densities to nodes
    for n,d in K.nodes_iter(data=True):
        d[density_attr] += d[SMOOTH_TMP]

def assemble_transcript_graph(G, strand, partial_paths, kmax, fraction_major_path, max_paths):
    """
    partial_path_data: dictionary where key is partial path tuple and value
    is density
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
    # kmax should not be set to less than '2'
    kmax = max(2, kmax)
    # sort partial path data by increasing length of the path
    sorted_partial_paths = sorted(partial_paths, key=lambda x: len(x[0])) 
    # choose 'k' no larger than largest partial path
    k = min(kmax, len(sorted_partial_paths[-1][0]))
    # append/prepend dummy nodes to start/end nodes
    start_nodes, end_nodes = add_dummy_start_end_nodes(G, k)
    # initialize k-mer graph
    K = create_kmer_graph(G, k)
    # add partial paths k-mers to graph
    add_partial_path_kmers(G, K, sorted_partial_paths, k)
    # smooth kmer graph
    smooth_kmer_graph(K)
    # find up to 'max_paths' paths through graph
    path_info_list = []
    for kmer_path, density in find_suboptimal_paths(K, fraction_major_path,
                                                    max_paths):
        # get nodes from kmer path
        path = list(kmer_path[0])
        path.extend(kmer[-1] for kmer in kmer_path[1:])
        path = [n for n in path if n.start >= 0]
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
    return path_info_list


