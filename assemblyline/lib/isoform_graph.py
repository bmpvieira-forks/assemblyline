'''
Created on Nov 23, 2010

@author: mkiyer
'''
import logging
import networkx as nx
import collections

from bx.intervals.intersection import Interval, IntervalTree
from bx.intervals.cluster import ClusterTree

from cNode import Node, strand_int_to_str, strand_str_to_int, node_type_to_str, cmp_strand, interval_overlap

EXON = 0
INTRON = 1
POS_STRAND = 0
NEG_STRAND = 1
NO_STRAND = 2

# a structure for managing isoform graph information
# stores the node plus its immediate edges and scores
NodeDataTuple = collections.namedtuple('NodeDataTuple', ['node', 'ids', 'predecessors', 'successors'])
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

def cluster_overlapping_nodes(nodes):
    # cluster nodes
    cluster_tree = ClusterTree(0,1)
    for i,node in enumerate(nodes):
        cluster_tree.insert(node.start, node.end, i)
    # find groups of overlapping nodes
    for start, end, indexes in cluster_tree.getregions():        
        yield start, end, [nodes[i] for i in indexes]

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
    if succs is not None:
        succs.update(orig_succs)
    for orig_node in orig_nodes:
        ids.update(G.node[orig_node]['ids'])
        preds.update(G.predecessors(orig_node)) 
        succs.update(G.successors(orig_node))
    return MergeTuple(newnode, frozenset(ids), frozenset(preds), frozenset(succs))

def get_overhang(leaf_node, match_node):
    start_diff = match_node.start - leaf_node.start
    end_diff = leaf_node.end - match_node.end
    return start_diff, end_diff

def get_free_ends(G, node):
    """returns whether the (left,right) genomic ends of the node are free ends without edges"""
    fiveprime = G.in_degree(node) == 0
    threeprime = G.out_degree(node) == 0
    if node.strand == NEG_STRAND:
        return threeprime, fiveprime
    else:
        return fiveprime, threeprime

def trim(G, leaf_node, match_node, overhang_threshold):
    left_overhang, right_overhang = get_overhang(leaf_node, match_node)
    leaf_left_free, leaf_right_free = get_free_ends(G, leaf_node)
    match_left_free, match_right_free = get_free_ends(G, match_node)
    # to perform trimming the candidate leaf must overhang into an intron-incompatible
    # area of the target transcript within a specified threshold amount
    trim_left = leaf_left_free and (not match_left_free) and left_overhang > 0 and left_overhang < overhang_threshold 
    trim_right = leaf_right_free and (not match_right_free) and right_overhang > 0 and right_overhang < overhang_threshold    
    logging.debug("trim left:free=%s,overhang=%d,trim=%s right:free=%s,overhang=%d,trim=%s" % 
                  (leaf_left_free, left_overhang, trim_left, leaf_right_free, right_overhang, trim_right))
    if trim_left or trim_right:
        start = match_node.start if trim_left else leaf_node.start
        end = match_node.end if trim_right else leaf_node.end
        newnode = Node(start, end, merge_strand(leaf_node.strand, match_node.strand), EXON)
        return newnode    
    else:
        return None 

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
    left_free, right_free = get_free_ends(G, node)
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
                    logging.debug("Found trimmed intron-exon join: %s %s" % (intron, node))
                    left_intron_dict[intron.end].add(intron)                    
        if right_free:
            for intron in intron_starts[node.end][strand]:
                right_intron_dict[node.end].add(intron)
            #right_intron_dict[node.end].update(intron_starts[node.end][strand])
            if do_trim:
                for interval in intron_start_tree[strand].find(node.end, node.end):
                    intron = interval.value 
                    logging.debug("Found trimmed intron-exon join: %s %s" % (intron, node))
                    right_intron_dict[intron.start].add(intron)
        if node.start not in left_intron_dict:
            left_intron_dict[node.start] = set()
        if node.end not in right_intron_dict:
            right_intron_dict[node.end] = set()
        for start, left_introns in left_intron_dict.iteritems():
            logging.debug("start=%d left_introns=%s" % (start, str(left_introns)))                
            for end, right_introns in right_intron_dict.iteritems():
                logging.debug("start=%d end=%d left_introns=%s right_introns=%s" % (start, end, str(left_introns), str(right_introns)))                
                if len(left_introns) == 0 and len(right_introns) == 0:
                    continue
                newnode = Node(start, end, strand, EXON)
                preds, succs = order_by_strand(left_introns, right_introns, strand)
                logging.debug("Found valid join %s -> %s -> %s" % (preds, newnode, succs))                
                yield newnode, node, preds, succs

def merge_node_pair(G, n1_exon, n1_introns, n2_exon, n2_introns, 
                    do_trim=False, overhang_threshold=0):
    n1_nodes = set([n1_exon])
    n2_nodes = set([n2_exon])
    if do_trim:
        # get trimmed node of n1 with respect to n2
        trimmed_node = trim(G, n1_exon, n2_exon, overhang_threshold)
        if trimmed_node is not None:
            n1_nodes.add(trimmed_node)
        # get trimmed node of n2 with respect to n1
        trimmed_node = trim(G, n2_exon, n1_exon, overhang_threshold)
        if trimmed_node is not None:
            n2_nodes.add(trimmed_node)
    for n1 in n1_nodes:
        for n2 in n2_nodes:
            # ensure that there are no introns interfering between the
            # two nodes
            if (is_intron_compatible(n1, n1_introns, n2, n2_introns)):
                logging.debug("Found intron compatible nodes:")
                logging.debug("NODE1:%s" % str(n1))
                logging.debug("NODE2:%s" % str(n2))
                start = n1.start if n1.start <= n2.start else n2.start
                end = n1.end if n1.end >= n2.end else n2.end                
                newnode = Node(start, end, merge_strand(n1.strand, n2.strand), EXON)
                yield newnode, [n1_exon, n2_exon]

def merge_nodes(G, leaf_nodes, internal_nodes, 
                do_trim=False, overhang_threshold=0):
    # build an interval tree to find overlapping node pairs
    nodes = list(leaf_nodes) + list(internal_nodes)
    exon_tree = IntervalTree()
    for i,node in enumerate(nodes):    
        exon_tree.insert_interval(Interval(node.start, node.end, value=i))
    # find nodes overlapping each leaf node
    for i,n1 in enumerate(leaf_nodes):
        n1_introns = G.predecessors(n1) + G.successors(n1)        
        for hit in exon_tree.find(n1.start-1, n1.end+1):
            # do not need to compare node against itself or any leaf nodes
            # with lower indices because merging A,B and B,A would be 
            # redundant and unnecessary
            j = hit.value
            if i >= j:
                continue
            n2 = nodes[j]
            n2_introns = G.predecessors(n2) + G.successors(n2)            
            for merged_node, original_nodes in merge_node_pair(G, n1, n1_introns, n2, n2_introns, do_trim, overhang_threshold):                
                yield merged_node, original_nodes


def cluster_single_exon_transcripts(G, trim, overhang_threshold):
    '''
    optimize the merging by first clustering single exon transcripts
    that do not overlap any multi-exon transcripts
    '''
    # build interval trees and find single exon
    # nodes that are candidates for clustering
    intron_tree = IntervalTree()
    exon_tree = IntervalTree()
    single_exon_nodes = []
    for node in G.nodes_iter():
        if node.node_type == INTRON:
            intron_tree.insert_interval(Interval(node.start, node.end, strand=node.strand))
        elif (node.strand != NO_STRAND) or (G.degree(node) > 0):
            exon_tree.insert_interval(Interval(node.start, node.end, strand=node.strand))            
        else:
            # TODO: NO_STRAND is not required here, but single exon nodes
            # that have strand information need to be separated and clustered
            # separately using separate "free_nodes" lists
            single_exon_nodes.append(node)
    intronlist_exon_map = collections.defaultdict(lambda: [])
    for node in single_exon_nodes:
        # should only cluster nodes that are not candidates for trimming
        # and eventual joining to graph
        if trim:
            overlapping_exons = exon_tree.find(node.start, node.end)
            can_trim_left = any((0 <= (e.start - node.start) <= overhang_threshold) for e in overlapping_exons)
            can_trim_right = any((0 <= (node.end - e.end) <= overhang_threshold) for e in overlapping_exons)
            if can_trim_left or can_trim_right:
                #logging.debug("pre-clustering skipping trim-able node=%s" % node)
                continue
        # organize exons by how they overlap introns by 
        # building an (intron list) -> (single exon node) mapping
        overlapping_introns = intron_tree.find(node.start, node.end)
        intron_tuple = tuple((intron.start, intron.end, intron.strand)
                             for intron in overlapping_introns)
        intronlist_exon_map[intron_tuple].append(node)
    # single exon transcript nodes that overlap identical 
    # intron lists can be directly clustered without sacrificing 
    # potential merging with multi-exon transcripts
    for intronlist, exonlist in intronlist_exon_map.iteritems():
        cluster_tree = ClusterTree(0, 1)
        for index, node in enumerate(exonlist):
            cluster_tree.insert(node.start, node.end, index)            
        for start, end, exonlist_indexes in cluster_tree.getregions():
            if len(exonlist_indexes) <= 1:
                continue
            # create a new node
            orig_nodes = set(exonlist[index] for index in exonlist_indexes)
            ids = set()
            for orig_node in orig_nodes:
                ids.update(G.node[orig_node]['ids'])
            # add the new clustered node to the graph
            new_node = Node(start, end, NO_STRAND, EXON)
            merge_tuple = MergeTuple(new_node, frozenset(ids), frozenset(), frozenset())
            logging.debug("pre-clustering node=%s original nodes=%s" % (new_node, orig_nodes))
            yield merge_tuple, orig_nodes

def thread_transcript(G, exon_tree, strand, exons, exon_index=0, path=None):
    if path == None:
        path = []
    if exon_index == len(exons):
        yield path
    else:
        e = exons[exon_index]     
        # intersect with exons to find potential path nodes
        # this is approximate and will be constrained further
        nodes = [hit.value for hit in exon_tree.find(e.start, e.end)
                 if (cmp_strand(hit.value.strand, strand) and
                     (e.start >= hit.value.start) and
                     (e.end <= hit.value.end))]
        valid_nodes = []
        for n in nodes:
            # make sure there is a valid path to node 'n'
            if len(path) > 0 and nx.shortest_path(G, path[-1], n) == False:
                continue
            # find which ends of the exon are flexible
            left_free = (((exon_index == 0) and (strand != NEG_STRAND)) or
                         ((exon_index == len(exons)) and (strand == NEG_STRAND)))
            right_free = (((exon_index == len(exons)) and (strand != NEG_STRAND)) or
                          ((exon_index == 0) and (strand == NEG_STRAND)))
            # ensure that exon-intron splice sites are constrained
            if (not left_free) and (e.start != n.start):
                continue
            if (not right_free) and (e.end != n.end):
                continue
            valid_nodes.append(n)
        if len(valid_nodes) == 0:
            # if there is no valid node to continue this path (which can 
            # happen due to trimming and other modifications to the original
            # transcripts) return this path early
            yield path
            # find new paths starting from the next node
            for newpath in thread_transcript(G, exon_tree, strand, exons, exon_index+1, path=None):
                yield newpath
        else:
            # found some valid paths
            for node in valid_nodes:  
                for newpath in thread_transcript(G, exon_tree, strand, exons, exon_index+1, path + [node]):
                    yield newpath

def get_connected_components(G):
    H = G.to_undirected()
    connected_nodes = nx.connected_components(H)
    for nodelist in connected_nodes:
        yield nx.DiGraph(G.subgraph(nodelist))

def find_paths_from_start(G, partial_path_dict, start_node, path=None):
    if path == None:
        path = []
    logging.debug('CURRENT_PATH=%s' % str(path))
    if G.out_degree(start_node) == 0:
        yield path + [start_node]
        logging.debug('END NODE REACHED=%s FINAL PATH=%s' % (str(start_node), str(path)))
    else:
        # search for partial paths starting from this node
        if start_node in partial_path_dict:
            for partial_path in partial_path_dict[start_node]:
                logging.debug('Adding partial path=%s' % str(partial_path))                
                for newpath in find_paths_from_start(G, partial_path_dict, partial_path[-1], path + list(partial_path[:-1])):
                    yield newpath
        else:
            logging.debug('ADDING NODE=%s' % str(start_node))
            # no partial path exists so use all possible paths
            for next_node in G.successors(start_node):
                logging.debug('NEXT NODE=%s' % str(next_node))
                for newpath in find_paths_from_start(G, partial_path_dict, next_node, path + [start_node]):
                    yield newpath


def find_paths(G, partial_paths):
    # make a dictionary so it is possible to look up 
    # partial paths based on the start node
    partial_path_dict = collections.defaultdict(lambda: [])
    for p in partial_paths:
        for i in xrange(len(p)-1):
            partial_path_dict[p[i]].append(p[i:])
    partial_path_dict = dict(partial_path_dict)        
    # break up the graph into connected subgraphs
    for Gsub in get_connected_components(G):
        # find start nodes
        start_nodes = [n for (n,d) in Gsub.in_degree_iter()
                       if (d == 0)]
        paths = []        
        for start_node in start_nodes:
            paths.extend(find_paths_from_start(Gsub, partial_path_dict, start_node))
        yield paths

def score_path(G, id_score_map, path):
    total_length = 0
    total_score = 0
    for node in path:
        total_length += (node.end - node.start)
        total_score += sum(id_score_map[id] for id in G.node[node]['ids'])
    return 1e3 * total_score / total_length

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
        logging.debug("CONNECTING EXONS")
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
        logging.debug("/CONNECTING EXONS")

    def _collapse_cluster(self, leaf_nodes, internal_nodes, 
                          merge_tuples, replaced_nodes,
                          trim=False, overhang_threshold=0):
        # find intron-compatible nodes that can be merged
        logging.debug("MERGING")
        #merged_node_dict = collections.defaultdict(lambda: set())
        merged_node_dict = collections.defaultdict(lambda: MergeInfo())
        for merged_node, orig_nodes in merge_nodes(self.G, leaf_nodes, internal_nodes, trim, overhang_threshold):
            merged_node_dict[merged_node].update_nodes(orig_nodes)         
#        for merged_node, orig_nodes in merged_node_dict.iteritems():
        for merged_node, merge_info in merged_node_dict.iteritems():
            replaced_nodes.update(merge_info.nodes)
            merge_tuple = make_merge_tuple(self.G, merged_node, merge_info.nodes, merge_info.preds, merge_info.succs)
            merge_tuples.append(merge_tuple)
            logging.debug("Merged intron compatible exons: %s" % str(merge_tuple))
        logging.debug("/MERGING")

    def collapse(self, trim=False, overhang_threshold=0):
        logging.info("COLLAPSE nodes=%d" % (len(self.G)))         
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
        # perform a pre-clustering step to collapse single
        # exon transcripts using a much simpler merge method
        for merge_tuple, replaced_nodes in cluster_single_exon_transcripts(self.G, trim, overhang_threshold):
            self.G.remove_nodes_from(replaced_nodes)
            self._add_merged_node(merge_tuple)
        logging.info("PRECLUSTERING nodes=%d" % (len(self.G)))
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
            logging.debug("cluster start=%d end=%d transcripts=%d" % (start, end, len(indexes)))
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
                    is_leaf = (self.G.in_degree(node) == 0 or self.G.out_degree(node) == 0)
                    if node.strand == NO_STRAND:
                        strands = [POS_STRAND, NEG_STRAND]
                    else:
                        strands = [node.strand]
                    for strand in strands:
                        if is_leaf:
                            strand_leaf_nodes[strand].append(node)
                        else:
                            strand_internal_nodes[strand].append(node)
                # collapse transcripts on each strand                
                merge_tuples = []
                replaced_nodes = set()            
                for strand in (POS_STRAND, NEG_STRAND):            
                    # get leaf nodes
                    leaf_nodes = strand_leaf_nodes[strand]
                    if len(leaf_nodes) == 0:
                        continue
                    # get internal nodes
                    internal_nodes = strand_internal_nodes[strand]
                    if len(leaf_nodes) + len(internal_nodes) < 2:
                        continue
                    logging.debug("cluster start=%d end=%d strand=%s leaf_nodes=%d internal_nodes=%d" % 
                                  (start, end, strand_int_to_str(strand), 
                                   len(leaf_nodes), len(internal_nodes)))
                    self._collapse_cluster(set(leaf_nodes), set(internal_nodes), 
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
                logging.debug("cluster start=%d end=%d iterations=%d merge operations performed=%d" % (start, end, iters, merges))
                # join free exon ends to compatible introns
                joins = 0
                merge_tuples = []
                replaced_nodes = set()
                for merge_tuple, orig_nodes in self._connect_leaf_exons(cluster_nodes,
                                                                        intron_starts, 
                                                                        intron_ends, 
                                                                        intron_start_tree, 
                                                                        intron_end_tree, 
                                                                        trim):
                    replaced_nodes.update(orig_nodes)
                    merge_tuples.append(merge_tuple)
                    logging.debug("JOIN tuple: %s" % str(merge_tuple))
                    logging.debug("REPLACED nodes: %s" % orig_nodes)
                    joins += 1                
                logging.debug("REPLACED nodes: %s" % replaced_nodes)
                self.G.remove_nodes_from(replaced_nodes)
                cluster_nodes.difference_update(replaced_nodes)
                logging.debug("ADDING nodes: %s" % str(merge_tuples))
                for merge_tuple in merge_tuples:
                    logging.debug("Adding JOIN: %s" % str(merge_tuple))
                    self._add_merged_node(merge_tuple)
                    cluster_nodes.add(merge_tuple.node)
                logging.debug("cluster start=%d end=%d iterations=%d join operations performed=%d" % (start, end, iters, joins))
                merges += joins
                iters += 1
        logging.debug("/COLLAPSE")
        
    def get_isoforms(self, transcripts):
        # build interval tree
        exon_tree = IntervalTree()
        for node in self.G.nodes_iter():
            if node.node_type == EXON:
                exon_tree.insert_interval(Interval(node.start, node.end, value=node))
        # find paths through graph that support initial transcripts
        logging.debug("THREADING")
        logging.debug("num transcripts=%d" % len(transcripts))        
        threaded_paths = set()
        id_score_map = {}
        for t in transcripts:
            # add scores to the score lookup table
            id_score_map[t.id] = t.score
            for e in t.exons:
                id_score_map[e.id] = e.score
            # ignore single exon transcripts
            if len(t.exons) == 1:
                continue
            # convert to integer strand
            strand = strand_str_to_int(t.strand)
            # flip exons for minus strand to keep path in correct order
            exons = t.exons[::-1] if strand == NEG_STRAND else t.exons
            for path in thread_transcript(self.G, exon_tree, strand, exons):
                # ignore paths of length 1 (useless for threading)
                if len(path) > 1:
                    # TODO: check path for now
                    for i in xrange(len(path)-1):
                        for j in xrange(i+1, len(path)):
                            if nx.shortest_path(self.G, path[i], path[j]) == False:
                                logging.error("PATH ERROR")
                                assert False
                    threaded_paths.add(tuple(path))
                for x in path:
                    assert x.node_type != INTRON
        logging.debug("threaded %s paths" % len(threaded_paths))
        for p in threaded_paths:
            logging.debug("threaded path=%s" % str([str(n) for n in p]))
            logging.debug("threaded node types=%s" % str([str(n.node_type) for n in p]))
            logging.debug("threaded exon sizes=%s" % str([str(n.end - n.start) for n in p]))
            logging.debug("threaded exon starts=%s" % str([str(n.start) for n in p]))
        logging.debug("/THREADING")
        # using transcript-threaded paths as a guide find full paths
        # through the graph        
        for gene_paths in find_paths(self.G, threaded_paths):
            # find non-redundant set of paths
            gene_paths = set(tuple(p) for p in gene_paths)
            scored_gene_paths = [(p, score_path(self.G, id_score_map, p)) for p in gene_paths]
            yield scored_gene_paths
