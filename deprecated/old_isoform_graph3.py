'''
Created on Dec 6, 2010

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

from cNode import Node, strand_int_to_str, strand_str_to_int, node_type_to_str, cmp_strand, interval_overlap
from base import EXON, INTRON, POS_STRAND, NEG_STRAND, NO_STRAND
from path_finder import thread_transcript, find_paths, score_path

# a structure for managing isoform graph information
# stores the node plus its immediate edges and scores
NodeEdgesTuple = collections.namedtuple('NodeEdgesTuple', ['node', 'predecessors', 'successors'])
MergeTuple = collections.namedtuple('MergeTuple', ['node', 'ids', 'predecessors', 'successors'])

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

def is_intron_compatible(exon1, introns1, exon2, introns2):
    if any(interval_overlap(exon1, intron) for intron in introns2):
        return False
    if any(interval_overlap(exon2, intron) for intron in introns1):
        return False
    return True
 
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

def get_overhang(leaf_node, match_node):
    start_diff = match_node.start - leaf_node.start
    end_diff = leaf_node.end - match_node.end
    return start_diff, end_diff

def get_free_ends(strand, in_degree, out_degree):
    """returns whether the (left,right) genomic ends of the node are free ends without edges"""
    fiveprime = (in_degree == 0)
    threeprime = (out_degree == 0)
    if strand == NEG_STRAND:
        return threeprime, fiveprime
    else:
        return fiveprime, threeprime

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

def separate_leaf_and_internal_nodes(G, nodes):
    # partition nodes by strand and whether node is a leaf or internal
    leaf_nodes = []
    internal_nodes = [] 
    for node in nodes:
        preds, succs = G.predecessors(node), G.successors(node)
        if len(preds) == 0 or len(succs) == 0:
            leaf_nodes.append(NodeEdgesTuple(node, preds, succs))
        else:
            internal_nodes.append(NodeEdgesTuple(node, preds, succs))
    return leaf_nodes, internal_nodes

def make_intron_indexes(G, do_trim=False, overhang_threshold=0):    
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
            if do_trim:
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

def trim(leaf_node, leaf_in_degree, leaf_out_degree,
         match_node, match_in_degree, match_out_degree,
         overhang_threshold):
    left_overhang, right_overhang = get_overhang(leaf_node, match_node)
    leaf_left_free, leaf_right_free = get_free_ends(leaf_node.strand, leaf_in_degree, leaf_out_degree)
    match_left_free, match_right_free = get_free_ends(match_node.strand, match_in_degree, match_out_degree)
    # to perform trimming the candidate leaf must overhang into an intron-incompatible
    # area of the target transcript within a specified threshold amount
    trim_left = leaf_left_free and (not match_left_free) and left_overhang > 0 and left_overhang < overhang_threshold 
    trim_right = leaf_right_free and (not match_right_free) and right_overhang > 0 and right_overhang < overhang_threshold    
    #logging.debug("trim left:free=%s,overhang=%d,trim=%s right:free=%s,overhang=%d,trim=%s" % 
    #              (leaf_left_free, left_overhang, trim_left, leaf_right_free, right_overhang, trim_right))
    if trim_left or trim_right:
        start = match_node.start if trim_left else leaf_node.start
        end = match_node.end if trim_right else leaf_node.end
        newnode = Node(start, end, merge_strand(leaf_node.strand, match_node.strand), EXON)
        return newnode    
    else:
        return None 

def merge_node_pair(G, n1_tuple, n2_tuple,
                    do_trim=False, overhang_threshold=0):
    '''
    n1_tuple: NodeEdgesTuple(node, predecessors, successors)
    n2_tuple: NodeEdgesTuple(node, predecessors, successors)
    do_trim: True to enable trimming
    overhang_threshold: max distance to allow trimming
    '''
    n1_nodes = set([n1_tuple.node])
    n2_nodes = set([n2_tuple.node])
    if do_trim:
        # get trimmed node of n1 with respect to n2
        trimmed_node = trim(n1_tuple.node, len(n1_tuple.predecessors), len(n1_tuple.successors),
                            n2_tuple.node, len(n2_tuple.predecessors), len(n2_tuple.successors),
                            overhang_threshold)
        if trimmed_node is not None:
            n1_nodes.add(trimmed_node)
        # get trimmed node of n2 with respect to n1
        trimmed_node = trim(n2_tuple.node, len(n2_tuple.predecessors), len(n2_tuple.successors),
                            n1_tuple.node, len(n1_tuple.predecessors), len(n1_tuple.successors),
                            overhang_threshold)
        if trimmed_node is not None:
            n2_nodes.add(trimmed_node)
    for n1 in n1_nodes:
        for n2 in n2_nodes:
            # ensure that there are no introns interfering between the
            # two nodes
            n1_introns = itertools.chain(n1_tuple.predecessors, n1_tuple.successors)
            n2_introns = itertools.chain(n2_tuple.predecessors, n2_tuple.successors)
            if (is_intron_compatible(n1, n1_introns, n2, n2_introns)):
                #logging.debug("Found intron compatible nodes:")
                #logging.debug("NODE1:%s" % str(n1))
                #logging.debug("NODE2:%s" % str(n2))
                start = n1.start if n1.start <= n2.start else n2.start
                end = n1.end if n1.end >= n2.end else n2.end                
                newnode = Node(start, end, merge_strand(n1.strand, n2.strand), EXON)
                yield newnode, [n1_tuple.node, n2_tuple.node]

def merge_nodes(G, leaf_node_tuples, internal_node_tuples, 
                do_trim=False, overhang_threshold=0):
    '''
    G: nx.DiGraph
    leaf_node_tuples: (node, predecessors, successors) tuples for each leaf
    internal_node_tuples: (node, predecessors, successors) tuples for each internal node
    do_trim: set to True to enabled node trimming to facilitate merging
    overhang_threshold: nodes separated by fewer than N bp can be trimmed
    '''
    # build an interval tree to find overlapping node pairs
    node_tuples = list(leaf_node_tuples) + list(internal_node_tuples)
    exon_tree = IntervalTree()
    for i,node_tuple in enumerate(node_tuples):
        node = node_tuple.node
        exon_tree.insert_interval(Interval(node.start, node.end, value=i))
    # find nodes overlapping each leaf node
    for i,n1_tuple in enumerate(leaf_node_tuples):
        for hit in exon_tree.find(n1_tuple.node.start-1, n1_tuple.node.end+1):
            # do not need to compare node against itself or any leaf nodes
            # with lower indices because merging A,B and B,A would be 
            # redundant and unnecessary
            j = hit.value
            if i >= j:
                continue
            n2_tuple = node_tuples[j]
            for merged_node, original_nodes in merge_node_pair(G, n1_tuple, n2_tuple,
                                                               do_trim, overhang_threshold):                
                yield merged_node, original_nodes

#def cluster_single_exon_transcripts(G, trim, overhang_threshold):
#    '''
#    optimize the merging by first clustering single exon transcripts
#    that do not overlap any multi-exon transcripts
#    '''
#    # build interval trees and find single exon
#    # nodes that are candidates for clustering
#    intron_tree = IntervalTree()
#    exon_tree = IntervalTree()
#    single_exon_nodes = []
#    for node in G.nodes_iter():
#        if node.node_type == INTRON:
#            intron_tree.insert_interval(Interval(node.start, node.end, strand=node.strand))
#        elif (node.strand != NO_STRAND) or (G.degree(node) > 0):
#            exon_tree.insert_interval(Interval(node.start, node.end, strand=node.strand))            
#        else:
#            # TODO: NO_STRAND is not required here, but single exon nodes
#            # that have strand information need to be separated and clustered
#            # separately using separate "free_nodes" lists
#            single_exon_nodes.append(node)
#    intronlist_exon_map = collections.defaultdict(lambda: [])
#    for node in single_exon_nodes:
#        # should only cluster nodes that are not candidates for trimming
#        # and eventual joining to graph
#        if trim:
#            overlapping_exons = exon_tree.find(node.start, node.end)
#            can_trim_left = any((0 <= (e.start - node.start) <= overhang_threshold) for e in overlapping_exons)
#            can_trim_right = any((0 <= (node.end - e.end) <= overhang_threshold) for e in overlapping_exons)
#            if can_trim_left or can_trim_right:
#                #logging.debug("pre-clustering skipping trim-able node=%s" % node)
#                continue
#        # organize exons by how they overlap introns by 
#        # building an (intron list) -> (single exon node) mapping
#        overlapping_introns = intron_tree.find(node.start, node.end)
#        intron_tuple = tuple((intron.start, intron.end, intron.strand)
#                             for intron in overlapping_introns)
#        intronlist_exon_map[intron_tuple].append(node)
#    # single exon transcript nodes that overlap identical 
#    # intron lists can be directly clustered without sacrificing 
#    # potential merging with multi-exon transcripts
#    for intronlist, exonlist in intronlist_exon_map.iteritems():
#        cluster_tree = ClusterTree(0, 1)
#        for index, node in enumerate(exonlist):
#            cluster_tree.insert(node.start, node.end, index)            
#        for start, end, exonlist_indexes in cluster_tree.getregions():
#            if len(exonlist_indexes) <= 1:
#                continue
#            # create a new node
#            orig_nodes = set(exonlist[index] for index in exonlist_indexes)
#            ids = set()
#            for orig_node in orig_nodes:
#                ids.update(G.node[orig_node]['ids'])
#            # add the new clustered node to the graph
#            new_node = Node(start, end, NO_STRAND, EXON)
#            merge_tuple = MergeTuple(new_node, frozenset(ids), frozenset(), frozenset())
#            logging.debug("pre-clustering node=%s original nodes=%s" % (new_node, orig_nodes))
#            yield merge_tuple, orig_nodes

def cluster_transcripts(G, trim, overhang_threshold):
    '''
    optimize the merging by first clustering transcripts within intron compatible "domains"
    '''
    # find all introns    
    introns = [n for n in G.nodes_iter() if n.node_type == INTRON]
    cluster_start = min(n.start for n in G.nodes_iter())
    cluster_end = max(n.end for n in G.nodes_iter())    
    # add all intron start/end positions to a set
    intron_boundaries = {POS_STRAND: set(), NEG_STRAND: set()}
    for intron in introns:
        assert intron.strand != NO_STRAND
        intron_boundaries[intron.strand].add(intron.start)
        intron_boundaries[intron.strand].add(intron.end)
    # sort the intron boundary positions and add them to interval trees
    intron_boundary_trees = {POS_STRAND: IntervalTree(),
                             NEG_STRAND: IntervalTree()}
    for strand, intron_boundary_set in intron_boundaries.iteritems():        
        sorted_intron_boundaries = [cluster_start] + sorted(intron_boundary_set) + [cluster_end]
        tree = intron_boundary_trees[strand]
        domain_id = 0
        for start, end in itertools.izip(sorted_intron_boundaries[:-1], sorted_intron_boundaries[1:]):
            tree.insert_interval(Interval(start, end, strand=strand, value=domain_id))
            domain_id += 1
    # place exons within intron boundary domains
    single_domain_exon_map = {POS_STRAND: collections.defaultdict(lambda: []),
                              NEG_STRAND: collections.defaultdict(lambda: [])}
    multi_domain_exon_map = {POS_STRAND: collections.defaultdict(lambda: []),
                             NEG_STRAND: collections.defaultdict(lambda: [])}
    for exon in G.nodes_iter():
        if exon.node_type != EXON:
            continue
        if exon.strand == NO_STRAND:
            strands = (POS_STRAND, NEG_STRAND)
        else:
            strands = (exon.strand,)
        for strand in strands:
            overlapping_domains = intron_boundary_trees[strand].find(exon.start, exon.end)
            assert len(overlapping_domains) > 0
            if len(overlapping_domains) == 1:
                domain_interval = overlapping_domains[0]
                assert (exon.start >= domain_interval.start) and (exon.end <= domain_interval.end)
                single_domain_exon_map[strand][domain_interval.value].append(exon)
            else:
                for domain_interval in overlapping_domains:
                    multi_domain_exon_map[strand][domain_interval.value].append(exon)
    merge_tuples = []
    replaced_nodes = set()
    for strand, domain_exon_map in single_domain_exon_map.iteritems():
        # cluster exons that lie within a single intron boundary domain
        for domain_id,exonlist in domain_exon_map.iteritems():
            # add to cluster tree
            cluster_tree = ClusterTree(0, 1)
            for index, exon in enumerate(exonlist):
                cluster_tree.insert(exon.start, exon.end, index)
            # get clustered regions
            cluster_merge_tuples = []
            cluster_replaced_nodes = set()
            for start, end, exonlist_indexes in cluster_tree.getregions():
                # no need to process singleton clusters
                if len(exonlist_indexes) <= 1:
                    continue
                # create a new merged node
                orig_nodes = set(exonlist[index] for index in exonlist_indexes)
                newstrand = reduce(merge_strand, [n.strand for n in orig_nodes])                
                newnode = Node(start, end, newstrand, EXON)
                cluster_merge_tuples.append(make_merge_tuple(G, newnode, orig_nodes))
                cluster_replaced_nodes.update(orig_nodes)
                logging.debug("intron-domain start=%d end=%d strand=%d nodes=%d" % (start, end, strand, len(exonlist_indexes)))
            # before replacing nodes in this domain, must ensure they cannot merge
            # within nodes that overlap multiple domains
            multi_exonlist = multi_domain_exon_map[strand][domain_id]
            multi_leaf_nodes, multi_internal_nodes = separate_leaf_and_internal_nodes(G, multi_exonlist)
            single_leaf_nodes, single_internal_nodes = separate_leaf_and_internal_nodes(G, cluster_replaced_nodes)
            node_sets_to_merge = [(multi_leaf_nodes, single_leaf_nodes), 
                                  (multi_leaf_nodes, single_internal_nodes),
                                  (single_leaf_nodes, multi_internal_nodes)]

            # find intron-compatible nodes that can be merged
            merged_node_dict = collections.defaultdict(lambda: set())
            for set1_node_tuples, set2_node_tuples in node_sets_to_merge:
                #for n1 in set1_node_tuples:
                #    for n2 in set2_node_tuples:
                #        logging.debug("need to check if node %s can merge with node %s" % (n1, n2))
                for merged_node, orig_nodes in merge_nodes(G, set1_node_tuples, set2_node_tuples, trim, overhang_threshold):
                    merged_node_dict[merged_node].update(orig_nodes)
            for merged_node, orig_nodes in merged_node_dict.iteritems():
                replaced_nodes.update(orig_nodes)
                merge_tuple = make_merge_tuple(G, merged_node, orig_nodes)
                merge_tuples.append(merge_tuple)
            merge_tuples.extend(cluster_merge_tuples)
            replaced_nodes.update(cluster_replaced_nodes)
    return merge_tuples, replaced_nodes

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
        for neighbor in merge_tuple.predecessors:
            self.G.add_edge(neighbor, node)
        for neighbor in merge_tuple.successors:
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

    def _connect_leaf_exons(self, leaf_nodes, intron_starts, intron_ends, intron_start_tree, intron_end_tree, do_trim):
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
            for joined_node, orig_node, preds, succs in find_intron_exon_connections(self.G, node, intron_starts, intron_ends, 
                                                                                     intron_start_tree, intron_end_tree,
                                                                                     do_trim):
                join_info = join_info_dict[joined_node]
                join_info.nodes.add(orig_node)
                join_info.preds.update(preds)
                join_info.succs.update(succs)
        for joined_node, merge_info in join_info_dict.iteritems():
            merge_tuple = make_merge_tuple(self.G, joined_node, merge_info.nodes, merge_info.preds, merge_info.succs)
            yield merge_tuple, merge_info.nodes
        #logging.debug("/CONNECTING EXONS")

    def _collapse_cluster(self, leaf_node_tuples, internal_node_tuples, 
                          merge_tuples, replaced_nodes,
                          trim=False, overhang_threshold=0):
        # find intron-compatible nodes that can be merged
        #logging.debug("MERGING")
        #merged_node_dict = collections.defaultdict(lambda: set())
        merged_node_dict = collections.defaultdict(lambda: MergeInfo())
        for merged_node, orig_nodes in merge_nodes(self.G, leaf_node_tuples, internal_node_tuples, trim, overhang_threshold):
            merged_node_dict[merged_node].update_nodes(orig_nodes)         
#        for merged_node, orig_nodes in merged_node_dict.iteritems():
        for merged_node, merge_info in merged_node_dict.iteritems():
            replaced_nodes.update(merge_info.nodes)
            merge_tuple = make_merge_tuple(self.G, merged_node, merge_info.nodes, merge_info.preds, merge_info.succs)
            merge_tuples.append(merge_tuple)
            #logging.debug("Merged intron compatible exons: %s" % str(merge_tuple))
        #logging.debug("/MERGING")

    def collapse(self, trim=False, overhang_threshold=0):
        logging.info("COLLAPSE called with nodes=%d" % (len(self.G)))
        # build intron intersection tables to add edges between
        # introns and exons
        intron_starts, intron_ends, intron_start_tree, intron_end_tree = make_intron_indexes(self.G, trim, overhang_threshold) 
        # now that all transcripts have been added, perform and initial
        # joining of free exon ends to compatible introns
        for merge_tuple, orig_nodes in self._connect_leaf_exons(self.G.nodes_iter(),
                                                                intron_starts, 
                                                                intron_ends, 
                                                                intron_start_tree, 
                                                                intron_end_tree, 
                                                                trim):
            # don't remove original nodes in this first step
            # because the nodes need the opportunity to merge
            # with all the rest of the nodes in the graph first
            self._add_merged_node(merge_tuple)
        # perform a pre-clustering step to collapse transcripts
        # within intron-compatible regions
        merge_tuples, replaced_nodes = cluster_transcripts(self.G, trim, overhang_threshold)
        self.G.remove_nodes_from(replaced_nodes)
        for merge_tuple in merge_tuples:
            self._add_merged_node(merge_tuple)
        logging.info("COLLAPSE clustered intron domain transcripts nodes=%d" % (len(self.G)))
        # perform a pre-clustering step to collapse single
        # exon transcripts using a much simpler merge method
        #for merge_tuple, replaced_nodes in cluster_single_exon_transcripts(self.G, trim, overhang_threshold):
        #    self.G.remove_nodes_from(replaced_nodes)
        #    self._add_merged_node(merge_tuple)
        # cluster exons
        cluster_tree = ClusterTree(0,1)
        exon_nodes = []
        exon_index = 0
        for node in self.G.nodes_iter():
            if node.node_type == EXON:
                # add to cluster tree
                cluster_tree.insert(node.start, node.end, exon_index)
                exon_nodes.append(node)
                exon_index += 1
        # find clusters of overlapping nodes
        for start, end, indexes in cluster_tree.getregions():
            # collapse each cluster of exon nodes independently
            logging.debug("COLLAPSE cluster start=%d end=%d nodes=%d" % (start, end, len(indexes)))
            # define initial starting set of cluster nodes
            cluster_nodes = set(exon_nodes[i] for i in indexes)
            merges = 1
            iters = 1
            while merges > 0:        
                # partition nodes by strand and whether node is a leaf or internal
                strand_leaf_nodes = {POS_STRAND: [], 
                                     NEG_STRAND: []}
                strand_internal_nodes = {POS_STRAND: [], 
                                         NEG_STRAND: []}
                for node in cluster_nodes:
                    preds, succs = self.G.predecessors(node), self.G.successors(node)
                    if node.strand == NO_STRAND:
                        strands = [POS_STRAND, NEG_STRAND]
                    else:
                        strands = [node.strand]
                    for strand in strands:
                        if len(preds) == 0 or len(succs) == 0:
                            strand_leaf_nodes[strand].append(NodeEdgesTuple(node, preds, succs))
                        else:
                            strand_internal_nodes[strand].append(NodeEdgesTuple(node, preds, succs))
                # collapse transcripts on each strand                
                merge_tuples = []
                replaced_nodes = set()            
                for strand in (POS_STRAND, NEG_STRAND):            
                    # get leaf nodes
                    leaf_node_tuples = strand_leaf_nodes[strand]
                    if len(leaf_node_tuples) == 0:
                        continue
                    # get internal nodes
                    internal_node_tuples = strand_internal_nodes[strand]
                    if len(leaf_node_tuples) + len(internal_node_tuples) < 2:
                        continue
                    logging.debug("COLLAPSE cluster start=%d end=%d strand=%s leaf_nodes=%d internal_nodes=%d" % 
                                  (start, end, strand_int_to_str(strand), 
                                   len(leaf_node_tuples), len(internal_node_tuples)))
                    self._collapse_cluster(leaf_node_tuples, internal_node_tuples, 
                                           merge_tuples, replaced_nodes, 
                                           trim, overhang_threshold)
                # remove all exons that were replaced by merges
                self.G.remove_nodes_from(replaced_nodes)
                cluster_nodes.difference_update(replaced_nodes)
                # reconstruct the graph by adding all the merged exons
                for merge_tuple in merge_tuples:
                    self._add_merged_node(merge_tuple)
                    cluster_nodes.add(merge_tuple.node)
                merges = len(merge_tuples)
                logging.debug("COLLAPSE cluster start=%d end=%d iterations=%d merges=%d nodes=%d" % (start, end, iters, merges, len(cluster_nodes)))
                # join free exon ends to compatible introns
                joins = 0
                merge_tuples = []
                replaced_nodes = set()
                leaf_nodes = [n for n,d in self.G.in_degree_iter(cluster_nodes) if d == 0]
                leaf_nodes.extend([n for n,d in self.G.out_degree_iter(cluster_nodes) if d == 0])
                for merge_tuple, orig_nodes in self._connect_leaf_exons(leaf_nodes,
                                                                        intron_starts, 
                                                                        intron_ends, 
                                                                        intron_start_tree, 
                                                                        intron_end_tree, 
                                                                        trim):
                    replaced_nodes.update(orig_nodes)
                    merge_tuples.append(merge_tuple)
                    joins += 1                
                #logging.debug("REPLACED nodes: %s" % replaced_nodes)
                self.G.remove_nodes_from(replaced_nodes)
                cluster_nodes.difference_update(replaced_nodes)
                #logging.debug("ADDING nodes: %s" % str(merge_tuples))
                for merge_tuple in merge_tuples:
                    self._add_merged_node(merge_tuple)
                    cluster_nodes.add(merge_tuple.node)
                logging.debug("COLLAPSE cluster start=%d end=%d iterations=%d joins=%d nodes=%d" % (start, end, iters, joins, len(cluster_nodes)))
                merges += joins
                iters += 1
        logging.debug("/COLLAPSE")

