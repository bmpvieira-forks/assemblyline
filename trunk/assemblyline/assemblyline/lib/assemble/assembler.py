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

def get_path_from_kmers(kmers):
    path = list(kmers[0])

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
        return
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
    else:
        G.add_edge(dummy_nodes[-1], n)

def add_dummy_start_end_nodes(G, k):
    for n in G.nodes_iter():
        if G.in_degree(n) == 0:
            add_dummy_nodes(G, n, k, REV)
        if G.out_degree(n) == 0:
            add_dummy_nodes(G, n, k, FWD)

def enumerate_paths_reverse(G, seed, length):
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

def enumerate_paths_forward(G, seed, length):
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
    # extend path in both directions
    L = (k - len(seed_path))
    rev_paths = enumerate_paths_reverse(G, seed_path[0], L, k)
    fwd_paths = enumerate_paths_forward(G, seed_path[-1], L, k)
    # enumerate kmers
    kmers = {}
    for rp in rev_paths:
        #print 'REV', rp
        for fp in fwd_paths:
            #print 'FWD', fp
            for pos in xrange(L+1):
                kmer = rp[pos:] + seed_path + fp[:pos]
                if pos == 0:
                    kmers[kmer] = SMOOTH_REV
                elif pos == L:
                    kmers[kmer] = SMOOTH_FWD
                else:
                    kmers[kmer] = None
    return kmers

def smooth_iteration(K, density_attr, smooth_attr):
    nodes = nx.topological_sort(K)
    for u in nodes:
        ud = K.node[u]
        smooth_density = ud[smooth_attr]
        succ = K.successors(u)
        total_nbr_density = sum(K.node[v][density_attr] for v in succ)
        for v in succ:
            vd = K.node[v]
            frac = vd[density_attr]/float(total_nbr_density)
            vd[SMOOTH_TMP] += frac * smooth_density
            vd[smooth_attr] += frac * smooth_density

def smooth_kmer_graph(K, density_attr=NODE_DENSITY):
    smooth_iteration(K, density_attr, SMOOTH_FWD)
    # reverse and smooth again
    K.reverse(copy=False)
    smooth_iteration(K, density_attr, SMOOTH_REV)
    K.reverse(copy=False)
    # apply densities to nodes
    for n,d in K.nodes_iter(data=True):
        d[density_attr] += d[SMOOTH_TMP]

def create_kmer_graph(G, partial_paths, k):
    """
    partial_paths: must be sorted by smallest -> largest path fragment
    """
    # private function to add kmers to graph and sum density
    def add_kmer(K, n, density):
        if n not in K:
            K.add_node(n, attr_dict={NODE_DENSITY: 0.0,
                                     SMOOTH_FWD: 0.0,
                                     SMOOTH_REV: 0.0,
                                     SMOOTH_TMP: 0.0})
        nd = K.node[n]
        nd[NODE_DENSITY] += density
    # create a graph of k-mers
    K = nx.DiGraph()
    # add kmers to graph
    while len(partial_paths) > 0:
        # get next largest partial path from list
        path, density = partial_paths[-1]
        # stop iteration when path lengths drop below 'k'
        if len(path) < k:
            break
        # create kmers and add to graph
        u = path[0:k]
        add_kmer(K, u, density)
        # the first kmer should be "smoothed"
        K.node[u][SMOOTH_REV] += density
        for i in xrange(1, len(path)-(k-1)):
            v = path[i:i+k]
            add_kmer(K, v, density)
            K.add_edge(u,v)
            u = v
        # the last kmer should be "smoothed"
        K.node[u][SMOOTH_FWD] += density
        # pop path from list
        partial_paths.pop()
    # add small partial paths by first extending 
    # them to be k-mers
    while len(partial_paths) > 0:
        # get next largest partial path from list
        path, density = partial_paths.pop()
        # find all kmers that could match this path
        kmer_dict = get_partial_path_kmers(G, path, k)
        total_density = sum(K.node[kmer][NODE_DENSITY] for kmer in kmer_dict
                            if kmer in K)
        if total_density == 0:
            # no k-mers with this path exist, so apply an averaged density 
            # to all possible matching kmers
            avg_density = density / float(len(kmer_dict))
            for kmer,smooth_dir in kmer_dict.iteritems():
                # add these k-mers since they don't exist
                add_kmer(K, kmer, avg_density)
                if smooth_dir is not None:
                    K.node[kmer][smooth_dir] += avg_density
        else:
            # apply density proportionately to all matching kmers 
            for kmer,smooth_dir in kmer_dict.iteritems():
                if kmer not in K:
                    continue
                nd = K.node[kmer]
                frac = nd[NODE_DENSITY] / float(total_density)
                nd[NODE_DENSITY] += (frac * density)
                if smooth_dir is not None:
                    nd[smooth_dir] += (frac * density)
    # perform "smoothing" of kmer graph
    smooth_kmer_graph(K)
    return K


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
    # sort partial path data by increasing length of the path
    sorted_partial_paths = sorted(partial_paths, key=lambda x: len(x[0])) 
    # choose 'k' no larger than largest partial path
    k = min(kmax, len(sorted_partial_paths[-1][0]))
    print "K = ", k
    # append/prepend dummy nodes to start/end nodes
    add_dummy_start_end_nodes(G, k)
    # build a k-mer graph using partial paths
    K = create_kmer_graph(G, sorted_partial_paths, k)
    # find up to 'max_paths' paths through graph
    path_info_list = []
    for kmer_path, density in find_suboptimal_paths(K, fraction_major_path,
                                                    max_paths):
        # get nodes from kmer path
        print kmer_path
        #path = list(kmer_path[0])
        #path.extend(kmer[-1] for kmer in kmer_path[1:])
        #print "PATH", path
        continue
        # get tss_id from path
        path = list(path)
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


