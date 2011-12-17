'''
Created on Dec 7, 2010

@author: mkiyer
'''
'''
Created on Nov 23, 2010

@author: mkiyer
'''
import logging
import networkx as nx
import collections
import itertools

from bx.intervals.intersection import Interval, IntervalTree
from bx.intervals.cluster import ClusterTree

from cNode import Node, strand_int_to_str, strand_str_to_int, node_type_to_str, cmp_strand, interval_overlap, interval_overlap_threshold
from base import EXON, INTRON, POS_STRAND, NEG_STRAND, NO_STRAND
from kclique import kcliques
from path_finder import thread_transcript, find_paths, score_path

# a structure for managing isoform graph information
# stores the node plus its immediate edges
ExonTuple = collections.namedtuple('ExonTuple', ['node', 'preds', 'succs', 'exon_domains', 'intron_domains'])
MergeTuple = collections.namedtuple('MergeTuple', ['node', 'ids', 'preds', 'succs'])

class MergeInfo(object):
    __slots__ = ('nodes', 'preds', 'succs')
    def __init__(self):
        self.nodes = set()
        self.preds = set()
        self.succs = set()
    def update_nodes(self, nodes):
        self.nodes.update(nodes)
    def update(self, other):
        self.nodes.update(other.nodes)
        self.preds.update(other.preds)
        self.succs.update(other.succs)

def merge_strand(strand1, strand2):
    if strand1 == strand2:
        return strand1
    elif strand1 == NO_STRAND:
        return strand2
    elif strand2 == NO_STRAND:
        return strand1
    else:
        logging.error("Incompatible strands")        
        assert False

def make_merge_tuple(G, newnode, orig_nodes, orig_preds=None, orig_succs=None):
    preds = set()
    succs = set()
    ids = set()
    if orig_preds is not None:
        preds.update(orig_preds)
    if orig_succs is not None:
        succs.update(orig_succs)
    for orig_node in orig_nodes:
        ids.update(G.node[orig_node]['ids'])
        preds.update(G.predecessors(orig_node)) 
        succs.update(G.successors(orig_node))
    return MergeTuple(newnode, frozenset(ids), frozenset(preds), frozenset(succs))

def get_free_ends(strand, in_degree, out_degree):
    """returns whether the (left,right) genomic ends of the node are free ends without edges"""
    fiveprime = (in_degree == 0)
    threeprime = (out_degree == 0)
    if strand == NEG_STRAND:
        return threeprime, fiveprime
    else:
        return fiveprime, threeprime

def make_intron_domain_trees(G):
    '''
    define genomic intervals representing unique intron 'domains'
    optimize the merging by first clustering transcripts within intron compatible "domains"
    '''
    # find all introns    
    introns = [n for n in G.nodes_iter() if n.node_type == INTRON]
    cluster_start = min(n.start for n in G.nodes_iter())
    cluster_end = max(n.end for n in G.nodes_iter())    
    # add all intron start/end positions to a set
    intron_boundaries = [set(), set()]
    for intron in introns:
        assert intron.strand != NO_STRAND
        intron_boundaries[intron.strand].add(intron.start)
        intron_boundaries[intron.strand].add(intron.end)
    # sort the intron boundary positions and add them to interval trees
    intron_domain_trees = [IntervalTree(), IntervalTree()]
    for strand, intron_boundary_set in enumerate(intron_boundaries):   
        sorted_intron_boundaries = [cluster_start]
        sorted_intron_boundaries.extend(sorted(intron_boundary_set))
        sorted_intron_boundaries.append(cluster_end)
        tree = intron_domain_trees[strand]
        domain_id = 0
        for start, end in itertools.izip(sorted_intron_boundaries[:-1], sorted_intron_boundaries[1:]):
            tree.insert_interval(Interval(start, end, strand=strand, value=domain_id))
            domain_id += 1
    return intron_domain_trees

def trim_interval(start, end, trim_left, trim_right, overhang_threshold):
    # TODO: move to C code
    trim_start, trim_end = start, end
    if trim_left:
        trim_start = start + overhang_threshold
        if trim_start > end:
            trim_start = end
    if trim_right:
        trim_end = end - overhang_threshold
        if trim_end < start:
            trim_end = start
    if trim_end < trim_start:
        return trim_end, trim_start
    return trim_start, trim_end

def get_node_intron_domains(exon, preds, succs, intron_domain_tree, overhang_threshold):
    # trim exons before computing exon domains
    if overhang_threshold > 0:
        left_free, right_free = get_free_ends(exon.strand, len(preds), len(succs))
        start, end, = trim_interval(exon.start, exon.end, 
                                    left_free, right_free,
                                    overhang_threshold)
    else:
        start, end = exon.start, exon.end
    # transform exon coordinates to intron domain regions
    exon_domains = set(interval.value for interval in intron_domain_tree.find(start, end))
    # convert introns into domain regions
    intron_domains = set()
    for intron in itertools.chain(preds, succs):
        intron_domains.update(interval.value for interval in intron_domain_tree.find(intron.start, intron.end))
    return exon_domains, intron_domains

def get_cluster_intron_domains(G, cluster_exons, intron_domain_trees, overhang_threshold=0):
    # partition nodes by strand
    strand_exons = [[], []]
    for node in cluster_exons:
        preds, succs = G.predecessors(node), G.successors(node)
        # add unstranded nodes to both strands
        if node.strand == NO_STRAND:
            assert len(preds) == 0
            assert len(succs) == 0
            strands = (POS_STRAND, NEG_STRAND)
        else:
            strands = (node.strand,)
        # find intron domains for each node        
        for strand in strands:
            edomains, idomains = get_node_intron_domains(node, preds, succs, 
                                                         intron_domain_trees[strand],
                                                         overhang_threshold)
            strand_exons[strand].append(ExonTuple(node, preds, succs, edomains, idomains))
    return strand_exons

def order_leaf_before_internal_exons(exon_tuples):
    # order the exon tuples so that the leaf
    # exons come before the internal exons
    ordered_exon_tuples = []
    internal_exon_tuples = []
    for e in exon_tuples:
        if len(e.preds) == 0 or len(e.succs) == 0:
            ordered_exon_tuples.append(e)
        else:
            internal_exon_tuples.append(e)
    num_leaf_exons = len(ordered_exon_tuples)
    ordered_exon_tuples.extend(internal_exon_tuples)
    return num_leaf_exons, ordered_exon_tuples

def find_intron_compatible_nodes(num_leaf_exons, exon_tuples):
    Gmerge = nx.Graph()
    for i in xrange(num_leaf_exons):
        e1_tuple = exon_tuples[i]
        for j in xrange(i+1, len(exon_tuples)):
            e2_tuple = exon_tuples[j]
            if (len(e1_tuple.exon_domains.intersection(e2_tuple.intron_domains)) == 0 and
                len(e2_tuple.exon_domains.intersection(e1_tuple.intron_domains)) == 0):
                Gmerge.add_edge(i,j)
#                interval_overlap_threshold(e1_tuple.node, e2_tuple.node, 1)):
                #print 'compatible', i, j, e1_tuple.node, e2_tuple.node
    return Gmerge

def make_exon_merge_tuple(G, cluster_start, cluster_end, strand, exon_tuples):
    preds = set()
    succs = set()
    ids = set()
    strand = NO_STRAND
    for exon_tuple in exon_tuples:
        preds.update(exon_tuple.preds)
        succs.update(exon_tuple.succs)
        node = exon_tuple.node
        ids.update(G.node[node]['ids'])
        strand = merge_strand(strand, node.strand)
    if strand == NEG_STRAND:
        if len(preds) == 0:
            end = cluster_end
            #end = max(e.node.end for e in exon_tuples)
        else:
            end = min(n.start for n in preds)
        if len(succs) == 0:
            start = cluster_start
            #start = min(e.node.start for e in exon_tuples)
        else:
            start = max(n.end for n in succs)
    else:
        if len(preds) == 0:
            start = cluster_start
            #start = min(e.node.start for e in exon_tuples)
        else:
            start = max(n.end for n in preds)
        if len(succs) == 0:
            end = cluster_end
            #end = max(e.node.end for e in exon_tuples)
        else:
            end = min(n.start for n in succs)
    return MergeTuple(Node(start, end, strand, EXON), ids, preds, succs)


def merge_exons(G, strand, exon_tuples, k):
    num_leaf_exons, exon_tuples = order_leaf_before_internal_exons(exon_tuples)
    # build a 'mergability' Graph where an edge between two nodes
    # indicates that they are intron compatible and can be merged
    Gmerge = find_intron_compatible_nodes(num_leaf_exons, exon_tuples)
    # find cliques of nodes to merge    
    logging.info("finding cliques k=%d" % k)
    cliques = kcliques(Gmerge, k)
    logging.info("found %d cliques" % len(cliques))    
    for clique_indexes in cliques:
        if len(clique_indexes) < 2:
            continue
        # cluster nodes
        cluster_tree = ClusterTree(0,1)
        for i in clique_indexes:
            node = exon_tuples[i].node
            cluster_tree.insert(node.start, node.end, i)
        # find groups of overlapping nodes
        for start, end, overlapping_indexes in cluster_tree.getregions():        
            if len(overlapping_indexes) > 1:
                # merge the exons into a new node
                clique_exon_tuples = [exon_tuples[i] for i in overlapping_indexes]
                merge_tuple = make_exon_merge_tuple(G, start, end, strand, clique_exon_tuples)
                logging.debug("Found clique size=%d" % (len(clique_exon_tuples)))
                for x,e in enumerate(clique_exon_tuples):
                    logging.debug("node%d: %s" % (x,e))
                logging.debug("Merged node=%s" % str(merge_tuple))
                yield merge_tuple, clique_exon_tuples


def make_intron_indexes(G, overhang_threshold=0):    
    # build an index of intron edges so that 
    # free exon ends can be joined with compatible introns
    intron_starts = collections.defaultdict(lambda: collections.defaultdict(lambda: set()))
    intron_ends = collections.defaultdict(lambda: collections.defaultdict(lambda: set()))
    # if trimming is enabled nodes will be trimmed in order
    # to be connected to existing introns
    intron_start_tree = {POS_STRAND: IntervalTree(),
                         NEG_STRAND: IntervalTree()}
    intron_end_tree = {POS_STRAND: IntervalTree(),
                       NEG_STRAND: IntervalTree()}    
    for node in G.nodes_iter():
        if node.node_type == INTRON:
            # maintain intron intersection structures
            assert node.node_type == INTRON
            # no unstranded introns allowed
            assert node.strand != NO_STRAND
            # keep track of positions where introns can be joined to exons
            intron_starts[node.start][node.strand].add(node)
            intron_ends[node.end][node.strand].add(node)
            # if trimming enabled build an intron index to find
            # trimmed exons that can be joined to introns
            if overhang_threshold > 0:
                intron_start_tree[node.strand].insert_interval(Interval(node.start, node.start + overhang_threshold, value=node))
                intron_end_tree[node.strand].insert_interval(Interval(node.end - overhang_threshold, node.end, value=node))
    return intron_starts, intron_ends, intron_start_tree, intron_end_tree

def find_intron_exon_connections(G, node, intron_starts, intron_ends, 
                                 intron_start_tree, intron_end_tree,
                                 do_trim):
    left_free, right_free = get_free_ends(node.strand, G.in_degree(node), G.out_degree(node))
    if not (left_free or right_free):
        return
    order_by_strand = lambda x,y,s: (y,x) if s == NEG_STRAND else (x,y)
    strands = [POS_STRAND, NEG_STRAND] if node.strand == NO_STRAND else [node.strand]
    for strand in strands:
        left_intron_dict = collections.defaultdict(lambda: set())
        right_intron_dict = collections.defaultdict(lambda: set())
        if left_free:
            for intron in intron_ends[node.start][strand]:
                left_intron_dict[node.start].add(intron)
            #left_intron_dict[node.start].update(intron_ends[node.start][strand])
            # if trimming enabled allow exons to overlap introns within threshold
            if do_trim:
                for interval in intron_end_tree[strand].find(node.start, node.start):
                    intron = interval.value 
                    #logging.debug("Found trimmed intron-exon join: %s %s" % (intron, node))
                    left_intron_dict[intron.end].add(intron)                    
        if right_free:
            for intron in intron_starts[node.end][strand]:
                right_intron_dict[node.end].add(intron)
            #right_intron_dict[node.end].update(intron_starts[node.end][strand])
            if do_trim:
                for interval in intron_start_tree[strand].find(node.end, node.end):
                    intron = interval.value 
                    #logging.debug("Found trimmed intron-exon join: %s %s" % (intron, node))
                    right_intron_dict[intron.start].add(intron)
        if node.start not in left_intron_dict:
            left_intron_dict[node.start] = set()
        if node.end not in right_intron_dict:
            right_intron_dict[node.end] = set()
        for start, left_introns in left_intron_dict.iteritems():
            #logging.debug("start=%d left_introns=%s" % (start, str(left_introns)))                
            for end, right_introns in right_intron_dict.iteritems():
                #logging.debug("start=%d end=%d left_introns=%s right_introns=%s" % (start, end, str(left_introns), str(right_introns)))                
                if len(left_introns) == 0 and len(right_introns) == 0:
                    continue
                newnode = Node(start, end, strand, EXON)
                preds, succs = order_by_strand(left_introns, right_introns, strand)
                #logging.debug("Found valid join %s -> %s -> %s" % (preds, newnode, succs))                
                yield newnode, node, preds, succs

def connect_leaf_exons(G, leaf_nodes, intron_starts, intron_ends, intron_start_tree, intron_end_tree, do_trim):
    '''
    join all free exon ends to compatible introns
    
    intron_starts: dict key=intron start position value=intron node
    intron_ends: dict key=intron end position value=intron node
    intron_start_tree: IntervalTree with intron start positions + overhang threshold for trimming
    intron_end_tree: IntervalTree with intron end positions - overhang threshold for trimming        
    '''
    #logging.debug("CONNECTING EXONS")
    join_info_dict = collections.defaultdict(lambda: MergeInfo()) 
    for node in leaf_nodes:
        for joined_node, orig_node, preds, succs in find_intron_exon_connections(G, node, intron_starts, intron_ends, 
                                                                                 intron_start_tree, intron_end_tree,
                                                                                 do_trim):
            join_info = join_info_dict[joined_node]
            join_info.nodes.add(orig_node)
            join_info.preds.update(preds)
            join_info.succs.update(succs)
    for joined_node, merge_info in join_info_dict.iteritems():
        merge_tuple = make_merge_tuple(G, joined_node, merge_info.nodes, merge_info.preds, merge_info.succs)
        yield merge_tuple, merge_info.nodes
    #logging.debug("/CONNECTING EXONS")


class IsoformGraph(object):
    def __init__(self):
        self.G = nx.DiGraph()

    @staticmethod
    def from_transcripts(transcripts):
        g = IsoformGraph()
        for t in transcripts:
            g.add_transcript(t)
        return g

    def _add_node(self, n, ids):
        if n not in self.G:  
            self.G.add_node(n, ids=set())
        self.G.node[n]['ids'].update(ids)

    def _add_merged_node(self, merge_tuple):
        # add node to graph
        node = merge_tuple.node
        self._add_node(node, merge_tuple.ids)
        # connect all neighbors to this new exon
        for neighbor in merge_tuple.preds:
            self.G.add_edge(neighbor, node)
        for neighbor in merge_tuple.succs:
            self.G.add_edge(node, neighbor)

    def add_transcript(self, transcript):
        exons = transcript.exons
        # convert from string strand notation ("+", "-") to integer (0, 1)
        strand = strand_str_to_int(transcript.strand)
        # add the first exon to initialize the loop
        # (all transcripts must have at least one exon)
        e1_node = Node(exons[0].start, exons[0].end, strand, EXON)
        self._add_node(e1_node, [exons[0].id])
        for e2 in exons[1:]:
            # add intron node            
            intron = Node(e1_node.end, e2.start, strand, INTRON)
            self._add_node(intron, [transcript.id])
            # add exon
            e2_node = Node(e2.start, e2.end, strand, EXON) 
            self._add_node(e2_node, [e2.id])
            # add edges from exon -> intron -> exon
            if strand == NEG_STRAND:
                self.G.add_edge(e2_node, intron)
                self.G.add_edge(intron, e1_node)
            else:
                self.G.add_edge(e1_node, intron)
                self.G.add_edge(intron, e2_node)
            # continue loop
            e1_node = e2_node

    def collapse(self, trim=False, overhang_threshold=0, max_cliques=3):
        logging.info("COLLAPSE called with nodes=%d" % (len(self.G)))
        if not trim:
            overhang_threshold = 0
        # use intron boundaries to define genomic regions,
        # of 'intron compatible domains' where merging can
        # be performed without checking for intron 
        # compatibility        
        intron_domain_trees = make_intron_domain_trees(self.G)
        # build intron intersection tables to add edges between
        # introns and exons
        intron_starts, intron_ends, intron_start_tree, intron_end_tree = make_intron_indexes(self.G, overhang_threshold)         
        # partition the exons into overlapping clusters
        # so that each cluster can be solved independently
        cluster_tree = ClusterTree(0,1)
        exon_nodes = []
        exon_index = 0
        for node in self.G.nodes_iter():
            if node.node_type == EXON:
                # add to cluster tree
                cluster_tree.insert(node.start, node.end, exon_index)
                exon_nodes.append(node)
                exon_index += 1
        # process clusters of overlapping exon nodes
        for start, end, indexes in cluster_tree.getregions():
            # collapse each cluster of exon nodes independently
            logging.debug("cluster start=%d end=%d nodes=%d" % (start, end, len(indexes)))
            # define initial starting set of cluster nodes
            cluster_exons = set(exon_nodes[i] for i in indexes)
            ops = 1
            iters = 1
            while ops > 0:
                merge_tuples = []
                replaced_nodes = set()
                # perform joining of free exon ends to compatible introns
                leaf_exons = [n for n,d in self.G.in_degree_iter(cluster_exons) if d == 0]
                leaf_exons.extend([n for n,d in self.G.out_degree_iter(cluster_exons) if d == 0])
                for merge_tuple, orig_nodes in connect_leaf_exons(self.G,
                                                                  leaf_exons,
                                                                  intron_starts, 
                                                                  intron_ends, 
                                                                  intron_start_tree, 
                                                                  intron_end_tree, 
                                                                  overhang_threshold > 0):
                    merge_tuples.append(merge_tuple)
                    replaced_nodes.update(orig_nodes)
                # get 'intron compatibility' domains associated with each exon
                strand_exon_tuples = get_cluster_intron_domains(self.G, cluster_exons, 
                                                                intron_domain_trees,
                                                                overhang_threshold)
                for strand, exon_tuples in enumerate(strand_exon_tuples):
                    for merge_tuple, orig_exon_tuples in merge_exons(self.G, strand, exon_tuples, k=max_cliques):
                        merge_tuples.append(merge_tuple)
                        replaced_nodes.update(e.node for e in orig_exon_tuples)
                # enact merges
                self.G.remove_nodes_from(replaced_nodes)
                cluster_exons.difference_update(replaced_nodes)
                for merge_tuple in merge_tuples:
                    self._add_merged_node(merge_tuple)
                    cluster_exons.add(merge_tuple.node)
                ops = len(merge_tuples)
                iters += 1
                #if iters == 2:
                #    break                
        logging.debug("/COLLAPSE")