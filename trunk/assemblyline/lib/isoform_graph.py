'''
Created on Nov 23, 2010

@author: mkiyer
'''
import logging
import networkx as nx
import collections

from bx.intervals.intersection import Interval, IntervalTree
from bx.intervals.cluster import ClusterTree

EXON = 0
INTRON = 1

def node_type_to_str(node_type):
    if node_type == EXON:
        return "EXON"
    if node_type == INTRON:
        return "INTRON"
    logging.error("Unknown node_type %d" % node_type)
    assert False

POS_STRAND = 0
NEG_STRAND = 1
NO_STRAND = 2

_strand_str_to_int = {"+": POS_STRAND,
                      "-": NEG_STRAND,
                      ".": NO_STRAND}
_strand_int_to_str = dict((v,k) for k,v in _strand_str_to_int.items())

def node_strand_str_to_int(strand):
    return _strand_str_to_int[strand]
def node_strand_int_to_str(strand):
    return _strand_int_to_str[strand]
def cmp_strand(a, b):
    if a == b:
        return True
    if (a == NO_STRAND) or (b == NO_STRAND):
        return True
    return False

class Node(object):
    __slots__ = ('start', 'end', 'strand', 'node_type')
    def __init__(self, start, end, strand, node_type):
        self.start = start
        self.end = end
        self.strand = strand
        self.node_type = node_type
    def __str__(self):
        return ("%s-%s(%s)[%s]" % (self.start, self.end, node_strand_int_to_str(self.strand), node_type_to_str(self.node_type)[0])) 
    def __repr__(self):
        return ("<%s(start=%d end=%d strand=%s node_type=%s>" %
                (self.__class__.__name__, self.start, self.end, node_strand_int_to_str(self.strand), node_type_to_str(self.node_type)))
    def __eq__(self, other):
        return ((self.start == other.start) and
                (self.end == other.end) and
                (self.strand == other.strand) and
                (self.node_type == other.node_type))
    def __ne__(self, other):
        return not self.__eq__(other)
    def __hash__(self):
        return (self.start << 18) | (self.end << 3) | (self.strand << 1) | (self.node_type)

    def cmp_strand(self, other):
        if self.strand == other.strand:
            return True
        elif (self.strand == NO_STRAND) or (other.strand == NO_STRAND):
            return True
        return False

    def is_overlapping(self, other):
        return (self.start < other.end) and (other.start < self.end) 

# information for merging transcripts
MergeTuple = collections.namedtuple('MergeTuple', ['node', 'scores', 'predecessors', 'successors'])

def interval_overlap(a, b):
    return (a.start < b.end) and (b.start < a.end)
    #return (a[0] < b[1]) and (b[0] < a[1])

def is_intron_compatible(exon1, introns1, exon2, introns2):
    for intron in introns2:
        if interval_overlap(exon1, intron):
            return False
    for intron in introns1:
        if interval_overlap(exon2, intron):
            return False
    return True

def find_overlapping_nodes(nodes):
    # cluster nodes
    cluster_tree = ClusterTree(0,1)
    for i,node in enumerate(nodes):
        cluster_tree.insert(node.start, node.end, i)
    # find groups of overlapping nodes
    for start, end, indexes in cluster_tree.getregions():        
        yield [nodes[i] for i in indexes]

def find_intron_compatible_nodes(G, leaf_nodes, internal_nodes, node_proxies):
    Gcompat = nx.Graph()
    nodes = list(leaf_nodes) + list(internal_nodes)
    # compare each leaf to rest of nodes
    for i in xrange(len(leaf_nodes)):
        n1 = nodes[i]                
        n1_introns = G.predecessors(n1) + G.successors(n1)        
        for j in xrange(i+1, len(nodes)):
            n2 = nodes[j]
            #logging.debug("N1=%s N2=%s" % (n1, n2))
            # ensure proxy nodes are not considered for merging
            if (n2 in node_proxies[n1]) or (n1 in node_proxies[n2]):
                continue
            #logging.debug("passed proxy N1=%s N2=%s" % (n1, n2))
            # ensure nodes weren't generated as trimmings of the same
            # original parent node
            if len(node_proxies[n1].intersection(node_proxies[n2])) > 0:
                continue
            n2_introns = G.predecessors(n2) + G.successors(n2)
            #logging.debug("passed proxy2 N1=%s N2=%s" % (n1, n2))
            # ensure that there are no introns interfering between the
            # two nodes, and that there is not already a path between
            # the nodes in the current graph
            if (is_intron_compatible(n1, n1_introns, n2, n2_introns) and
                (nx.shortest_path(G, n1, n2) == False) and
                (nx.shortest_path(G, n2, n1) == False)):
                logging.debug("Found intron compatible nodes:")
                logging.debug("NODE1:%s" % str(n1))
                logging.debug("NODE2:%s" % str(n2))                
                Gcompat.add_edge(i,j)
    # find all the maximal cliques in the graph, which equate to the
    # intron-compatible node groups
    for indexes in nx.find_cliques(Gcompat):
        if len(indexes) > 1:
            # return groups of nodes that overlap in genomic space for merging
            for overlapping_nodes in find_overlapping_nodes([nodes[i] for i in indexes]):
                if len(overlapping_nodes) > 1:
                    yield overlapping_nodes

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

def merge_intron_compatible_exons(G, nodes):    
    logging.debug("Nodes to merge: %s" % nodes)
    newnode = Node(min([n.start for n in nodes]),
                   max([n.end for n in nodes]),
                   reduce(merge_strand, [n.strand for n in nodes]),
                   EXON)
    newscores = []
    predecessors = set()
    successors = set()
    for n in nodes:
        newscores.extend(G.node[n]['scores'].items())
        predecessors.update(G.predecessors(n))
        successors.update(G.successors(n))
    return MergeTuple(newnode, tuple(newscores), frozenset(predecessors), frozenset(successors))

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
    logging.debug("trim overhang threshold=%d" % (overhang_threshold))
    logging.debug("left overhang=%d right overhang=%d" % (left_overhang, right_overhang))
    logging.debug("left end free=%s right end free=%s" % (leaf_left_free, leaf_right_free))
    logging.debug("trim left=%s right=%s" % (trim_left, trim_right))
    if trim_left or trim_right:
        start = match_node.start if trim_left else leaf_node.start
        end = match_node.end if trim_right else leaf_node.end
        newnode = Node(start, end, merge_strand(leaf_node.strand, match_node.strand), EXON)
        return newnode    
    else:
        return None 

def find_trimmable_nodes(G, leaf_nodes, internal_nodes, overhang_threshold):
    # must convert to ordered list to ensure leaf nodes 
    # before internal nodes
    nodes = list(leaf_nodes) + list(internal_nodes)
    # compare each leaf to rest of nodes
    for i in xrange(len(leaf_nodes)):
        n1 = nodes[i]        
        # compare each node against all other nodes in case
        # it could trim against more than one node at a time.
        # need to find all possible trimmings and choose a 
        # consensus otherwise strange exponential behavior occurs
        # in special cases during collapse
        trimmed_starts = set()
        trimmed_ends = set()
        logging.debug("leaf=%s" % str(n1))
        for j in xrange(len(nodes)):
            # do not compare against self
            if i == j:
                continue
            n2 = nodes[j]           
            # nodes must overlap to trim
            if not interval_overlap(n1, n2):
                continue
            logging.debug("match=%s" % str(n2))
            # see if leaf node can be trimmed to become 
            # intron compatible with another node
            newnode = trim(G, n1, n2, overhang_threshold)
            if newnode is not None:
                trimmed_starts.add(newnode.start)
                trimmed_ends.add(newnode.end)
        for s in trimmed_starts:
            for e in trimmed_ends:
                newnode = Node(s, e, n1.strand, EXON)
                logging.debug("trimmed node=%s" % str(newnode))
                yield newnode, n1      

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

def score_path(G, path):
    total_length = 0
    total_score = 0
    for node in path:
        total_length += (node.end - node.start)
        score_dict = G.node[node]['scores']
        total_score += sum(score_dict.values())/len(score_dict)
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

    def _add_node(self, n, scores):
        if n not in self.G:            
            self.G.add_node(n, scores=collections.defaultdict(lambda: 0.0))
        for label,score in scores:
            self.G.node[n]['scores'][label] += score

    def add_transcript(self, transcript):
        exons = transcript.exons
        # convert from string strand notation ("+", "-") to integer (0, 1)
        strand = node_strand_str_to_int(transcript.strand)
        # add the first exon to initialize the loop
        # (all transcripts must have at least one exon)
        e1_node = Node(exons[0].start, exons[0].end, strand, EXON)
        self._add_node(e1_node, [(transcript.label, exons[0].score)])
        for e2 in exons[1:]:
            # add intron node            
            intron = Node(e1_node.end, e2.start, strand, INTRON)
            self._add_node(intron, [(transcript.label,transcript.score)]) 
            # add exon
            e2_node = Node(e2.start, e2.end, strand, EXON) 
            self._add_node(e2_node, [(transcript.label,e2.score)])
            # add edges from exon -> intron -> exon
            if strand == NEG_STRAND:
                self.G.add_edge(e2_node, intron)
                self.G.add_edge(intron, e1_node)
            else:
                self.G.add_edge(e1_node, intron)
                self.G.add_edge(intron, e2_node)
            # continue loop
            e1_node = e2_node

    def _add_merged_node(self, merge_tuple):
        # add node to graph
        node = merge_tuple.node
        self._add_node(node, merge_tuple.scores)
        # connect all neighbors to this new exon
        for neighbor in merge_tuple.predecessors:
            self.G.add_edge(neighbor, node)
        for neighbor in merge_tuple.successors:
            self.G.add_edge(node, neighbor)
    
    def _collapse_cluster(self, leaf_nodes, internal_nodes, merge_tuples, replaced_nodes, trim=False, overhang_threshold=0):
        node_proxies = collections.defaultdict(lambda: set())
        # only perform trimming if enabled
        if trim:
            # find overhanging exons that can be trimmed from graph
            logging.debug("TRIMMING")
            trim_merge_tuples = set()
            trimmed_node_dict = collections.defaultdict(lambda: set())
            for trimmed_node, leaf_node in find_trimmable_nodes(self.G, leaf_nodes, internal_nodes, overhang_threshold):
                trimmed_node_dict[trimmed_node].add(leaf_node)
            # build new trimmed nodes to temporarily add to the graph
            for trimmed_node, orig_leaf_nodes in trimmed_node_dict.iteritems():
                # if the new node is already part of the graph, then
                # the trimmed node will be "absorbed" into the graph
                # during the merging step
                # check if the new node is already part of the graph
                if trimmed_node not in self.G:
                    # all the trimmed nodes must be explicity 
                    # removed unless they already exist or are
                    # used to facilitate a merge
                    replaced_nodes.add(trimmed_node)
                    # the trimmable nodes are by definition leaf nodes 
                    # during merging
                    leaf_nodes.add(trimmed_node)
                else:
                    # the original leaf nodes must be removed from the
                    # graph because the trimmed node will be integrated
                    # into the graph without the explicit need for merging
                    replaced_nodes.update(orig_leaf_nodes)
                # temporarily add the trimmed node to the graph
                # so that it can participate in merging
                scores = []
                preds = set()
                succs = set()
                for leaf_node in orig_leaf_nodes:
                    scores.extend(self.G.node[leaf_node]['scores'].items())
                    preds.update(self.G.predecessors(leaf_node)) 
                    succs.update(self.G.successors(leaf_node))                            
                trim_merge_tuple = MergeTuple(trimmed_node, tuple(scores), frozenset(preds), frozenset(succs))
                trim_merge_tuples.add(trim_merge_tuple)
                # maintain a link from the trimmed node to the original 
                # leaf nodes because all of them must be removed if any 
                # one of them participates in a merge
                node_proxies[trimmed_node].update(orig_leaf_nodes)
            # add the trimmed nodes to the graph
            for trim_merge_tuple in trim_merge_tuples:
                self._add_merged_node(trim_merge_tuple)
            logging.debug("/TRIMMING")
        # find and merge intron-compatible nodes
        logging.debug("MERGING")
        # find intron-compatible nodes that can be merged
        for nodes_to_merge in find_intron_compatible_nodes(self.G, leaf_nodes, internal_nodes, node_proxies):
            merge_tuple = merge_intron_compatible_exons(self.G, nodes_to_merge)
            logging.debug("Merged intron compatible exons: %s" % str(merge_tuple))
            merge_tuples.append(merge_tuple)
            replaced_nodes.update(nodes_to_merge)
            for n in nodes_to_merge:
                proxies = node_proxies[n]
                if len(proxies) > 0:
                    logging.debug("Removing proxy nodes: %s" % (proxies))            
                    replaced_nodes.update(proxies)
        logging.debug("/MERGING")

    def collapse(self, trim=False, overhang_threshold=0):
        # cluster exons
        cluster_tree = ClusterTree(0,1)
        nodes = []
        exon_index = 0
        for node in self.G.nodes_iter():
            if node.node_type == INTRON:
                continue
            # add to cluster tree
            cluster_tree.insert(node.start, node.end, exon_index)
            nodes.append(node)
            exon_index += 1
        # collapse each cluster of exon nodes independently
        logging.debug('COLLAPSE')
        # find clusters of overlapping nodes
        for start, end, indexes in cluster_tree.getregions():
            logging.debug("cluster start=%d end=%d transcripts=%d" % (start, end, len(indexes)))
            # define initial starting set of cluster nodes
            cluster_nodes = set(nodes[i] for i in indexes)
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
                    # get leaf nodes and internal nodes
                    leaf_nodes = strand_leaf_nodes[strand]
                    if len(leaf_nodes) == 0:
                        continue
                    internal_nodes = strand_internal_nodes[strand]
                    if len(leaf_nodes) + len(internal_nodes) < 2:
                        continue
                    logging.debug("cluster start=%d end=%d strand=%s leaf_nodes=%d internal_nodes=%d" % 
                                  (start, end, node_strand_int_to_str(strand), 
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
                iters += 1
        logging.debug("/COLLAPSE")

#    def _collapse_once(self, trim=False, overhang_threshold=0):
#        nodes = []
#        leaf_node_indexes = set()
#        # cluster exons by strand
#        cluster_trees = {POS_STRAND: ClusterTree(0,1), NEG_STRAND: ClusterTree(0,1)}
#        exon_index = 0
#        for node in self.G.nodes_iter():
#            if node.node_type == INTRON:
#                continue
#            if self.G.in_degree(node) == 0 or self.G.out_degree(node) == 0:
#                leaf_node_indexes.add(exon_index)
#            # add to cluster tree based on strand
#            if node.strand == NO_STRAND:
#                strands = [POS_STRAND, NEG_STRAND]
#            else:
#                strands = [node.strand]
#            for strand in strands:
#                cluster_trees[strand].insert(node.start, node.end, exon_index)
#            nodes.append(node)
#            exon_index += 1
#        # collapse each cluster of exon nodes independently
#        logging.debug('COLLAPSE')
#        merge_tuples = []
#        replaced_nodes = set()
#        for strand, cluster_tree in cluster_trees.iteritems():
#            logging.debug("strand=%s" % node_strand_int_to_str(strand))
#            # find clusters of overlapping nodes
#            for start, end, indexes in cluster_tree.getregions():
#                # get lists of leaf nodes and internal nodes
#                leaf_nodes = [nodes[i] for i in set(indexes).intersection(leaf_node_indexes)]                
#                logging.debug("leaf_nodes=%s" % leaf_nodes)
#                if len(leaf_nodes) == 0:
#                    continue
#                internal_nodes = [nodes[i] for i in set(indexes).difference(leaf_node_indexes)]
#                logging.debug("internal_nodes=%s" % internal_nodes)
#                if len(leaf_nodes) + len(internal_nodes) < 2:
#                    continue
#                self._collapse_cluster(set(leaf_nodes), set(internal_nodes), merge_tuples, replaced_nodes, trim, overhang_threshold)
#        logging.debug("/COLLAPSE")
#        # remove all exons that were replaced by merges
#        self.G.remove_nodes_from(replaced_nodes)
#        # reconstruct the graph by adding all the merged exons
#        for merge_tuple in merge_tuples:
#            self._add_merged_node(merge_tuple)
#        return len(merge_tuples)

#    def collapse(self, trim=False, overhang_threshold=0):        
#        merges = 1
#        iters = 1
#        while merges > 0:
#            #merges = self._collapse_once(trim=trim, overhang_threshold=overhang_threshold)
#            merges = self._collapse_once(trim=trim, overhang_threshold=overhang_threshold)
#            logging.debug("iterations=%d merge operations performed=%d" % (iters, merges))
#            iters += 1
#            #if iters == 2:
#            #    break

    def convert_introns_to_edges(self):
        H = nx.DiGraph()
        for n,d in self.G.nodes_iter(data=True):
            if n.node_type == INTRON:
                for pred in self.G.predecessors_iter(n):
                    for succ in self.G.successors(n):
                        H.add_node(pred, scores=self.G.node[pred]['scores'])
                        H.add_node(succ, scores=self.G.node[succ]['scores'])
                        H.add_edge(pred, succ, scores=d)
        
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
        for t in transcripts:
            # ignore single exon transcripts
            if len(t.exons) == 1:
                continue
            # convert to integer strand
            strand = node_strand_str_to_int(t.strand)
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
            scored_gene_paths = [(p, score_path(self.G, p)) for p in gene_paths]
            yield scored_gene_paths

#def find_intron_compatible_pair(G, n1, n1_introns, n2, n2_introns, 
#                                do_trim=False, overhang_threshold=0):
#    # first test the full length nodes
#    if is_intron_compatible(n1, n1_introns, n2, n2_introns):
#        return n1,n2
#    # search for trimmed matches if trimming is enabled
#    trimmed_nodes = []
#    if do_trim:
#        # can trim n1
#        n1_trimmed = trim(G, n1, n2, overhang_threshold)
#        if n1_trimmed is not None:
#            trimmed_nodes.append((n1_trimmed.end - n1_trimmed.start, n1_trimmed, n2))
#        # can trim n2
#        n2_trimmed = trim(G, n2, n1, overhang_threshold)
#        if n2_trimmed is not None:
#            trimmed_nodes.append((n2_trimmed.end - n2_trimmed.start, n2_trimmed, n1))
#        # select the largest trimmed node
#        if len(trimmed_nodes) > 0:
#            largest_trimmed_node = sorted(trimmed_nodes, key=operator.itemgetter(0))[0][1]
#            return largest_trimmed_node
#    return None
#                  
#                # find disconnected exons that can be connected
#                # to compatible introns
#                logging.debug("JOINING EDGES")
#                for node in leaf_nodes:
#                    predecessors = self.G.predecessors(node)        
#                    successors = self.G.successors(node)                        
#                    scores = self.G.node[node]['scores'].items(),
#                    neighbor_introns = set(self.G.neighbors(node)) 
#                    if node.strand == NO_STRAND:
#                        strands = [POS_STRAND, NEG_STRAND]
#                    else:
#                        strands = [node.strand]                    
#                    # look for a compatible intron edge
#                    for strand in strands:
#                        # check if this node begins at the
#                        # same point an intron ends
#                        for intron in intron_ends[(node.start, strand)].difference(neighbor_introns):
#                            if strand == NEG_STRAND:
#                                join_tuples.append((node, intron))
#                            else:
#                                join_tuples.append((intron, node))
#                            logging.debug("Joined exon=%s with intron=%s" % (node, intron))
#                        # check if this node ends at the same
#                        # point an intron begins 
#                        for intron in intron_starts[(node.end, strand)].difference(neighbor_introns):
#                            if strand == NEG_STRAND:
#                                join_tuples.append((intron, node))
#                            else:
#                                join_tuples.append((node, intron))
#                            logging.debug("Joined exon=%s with intron=%s" % (node, intron))
#                logging.debug("/JOINING EDGES")
#    def _collapse_once(self, trim=False, overhang_threshold=0):
#        # build interval tree
#        exon_tree = IntervalTree()        
#        for node in self.G.nodes_iter():
#            if node.node_type == EXON:
#                exon_tree.insert_interval(Interval(node.start, node.end, value=node))
#
##        # build an index of intron edges so that 
##        # free exon ends can be joined with compatible introns
##        intron_starts = collections.defaultdict(lambda: set())
##        intron_ends = collections.defaultdict(lambda: set())
##        for node in self.G.nodes_iter():
##            if node.node_type == INTRON:
##                intron_starts[(node.start, node.strand)].add(node)
##                intron_ends[(node.end, node.strand)].add(node)
#                
#        leaf_node_indexes = set()
#        nodes = []
#        # cluster exons by strand
#        cluster_trees = {POS_STRAND: ClusterTree(0,1), NEG_STRAND: ClusterTree(0,1)}
#        exon_index = 0
#        for node in self.G.nodes_iter():
#            if node.node_type == INTRON:
#                continue
#            if self.G.in_degree(node) == 0 or self.G.out_degree(node) == 0:
#                leaf_node_indexes.add(exon_index)
#            # add to cluster tree based on strand
#            if node.strand == NO_STRAND:
#                strands = [POS_STRAND, NEG_STRAND]
#            else:
#                strands = [node.strand]
#            for strand in strands:
#                cluster_trees[strand].insert(node.start, node.end, exon_index)
#            nodes.append(node)
#            exon_index += 1
#        
#        logging.debug('COLLAPSE')
#        merge_tuples = []
#        join_tuples = []
#        replaced_nodes = set()                
#        for strand, cluster_tree in cluster_trees.iteritems():
#            logging.debug("strand=%s" % node_strand_int_to_str(strand))
#            # find clusters of overlapping nodes
#            for start, end, indexes in cluster_tree.getregions():
#                # get lists of leaf nodes and internal nodes
#                leaf_nodes = [nodes[i] for i in set(indexes).intersection(leaf_node_indexes)]                
#                logging.debug("leaf_nodes=%s" % leaf_nodes)
#                if len(leaf_nodes) == 0:
#                    continue
#                internal_nodes = [nodes[i] for i in set(indexes).difference(leaf_node_indexes)]
#                logging.debug("internal_nodes=%s" % internal_nodes)
#                if len(leaf_nodes) + len(internal_nodes) < 2:
#                    continue
#                # only perform trimming if enabled
#                node_proxies = collections.defaultdict(lambda: set())
#                trimmed_nodes = set()
#                trim_merge_tuples = []
#                if trim:
#                    # find overhanging exons that can be trimmed from graph
#                    logging.debug("TRIMMING")
#                    for trimmed_node, leaf_node in find_trimmable_nodes(self.G, leaf_nodes, internal_nodes, overhang_threshold):
#                        # set unstranded trimmed node to the current
#                        # strand being collapsed before testing if 
#                        # node already exists in the graph
#                        if trimmed_node.strand == NO_STRAND:
#                            trimmed_node.strand = strand
#                        # if the new node is already part of the graph, then
#                        # essentially the trimmed node can be "absorbed" into
#                        # the graph without any additional merging
#                        # check if the new node is already part of the graph
#                        if trimmed_node in self.G:
#                            # the original leaf node can be removed from graph
#                            # because it can be trimmed to be equal to a different
#                            # node
#                            replaced_nodes.add(leaf_node)
#                        else:
#                            # temporarily add the trimmed node to the graph
#                            # so that it can participate in merging
#                            scores = self.G.node[leaf_node]['scores'].items()
#                            trim_merge_tuple = MergeTuple(trimmed_node, 
#                                                          tuple(scores),
#                                                          frozenset(self.G.predecessors(leaf_node)), 
#                                                          frozenset(self.G.successors(leaf_node)))
#                            trim_merge_tuples.append(trim_merge_tuple)
#                            # maintain a link from the trimmed node to the original 
#                            # leaf nodes so that all can be removed if any one of 
#                            # them participates in a merge
#                            node_proxies[trimmed_node].add(leaf_node)
#                            # the trimmable nodes are by definition leaf nodes 
#                            # during merging
#                            if trimmed_node not in trimmed_nodes:
#                                leaf_nodes.append(trimmed_node)
#                                trimmed_nodes.add(trimmed_node)
#                            # all the trimmed nodes will be explicity 
#                            # removed unless they are used to facilitate a merge
#                            replaced_nodes.add(trimmed_node)
#                    for trim_merge_tuple in trim_merge_tuples:
#                        # add the trimmable node to the graph
#                        self._add_merged_node(trim_merge_tuple)
#                    logging.debug("/TRIMMING")
#                # find and merge intron-compatible nodes
#                logging.debug("MERGING")
#                # find intron-compatible nodes that can be merged
#                for nodes_to_merge in find_intron_compatible_nodes(self.G, leaf_nodes, internal_nodes, node_proxies):
#                    merge_tuple = merge_intron_compatible_exons(self.G, nodes_to_merge)
#                    logging.debug("Merged intron compatible exons: %s" % str(merge_tuple))
#                    merge_tuples.append(merge_tuple)
#                    replaced_nodes.update(nodes_to_merge)
#                    for n in nodes_to_merge:
#                        logging.debug("Removing proxy nodes: %s" % (node_proxies[n]))            
#                        replaced_nodes.update(node_proxies[n])
#                logging.debug("/MERGING")              
#
#        logging.debug("/COLLAPSE")
#        # remove all exons that were replaced by merges
#        self.G.remove_nodes_from(replaced_nodes)
#        # reconstruct the graph by adding all the merged exons
#        for merge_tuple in merge_tuples:
#            self._add_merged_node(merge_tuple)
#        # add addition exon-intron edges
#        for from_node,to_node in join_tuples:
#            self.G.add_edge(from_node, to_node)
#        return len(merge_tuples) + len(join_tuples)

