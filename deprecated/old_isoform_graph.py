'''
Created on Nov 3, 2010

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

class Node(object):
    __slots__ = ('start', 'end', 'strand', 'node_type')
    def __init__(self, start, end, strand, node_type):
        self.start = start
        self.end = end
        self.strand = strand
        self.node_type = node_type
    def __str__(self):
        return ("%s-%s_%s_%s" % (self.start, self.end, self.strand, node_type_to_str(self.node_type))) 
    def __repr__(self):
        return ("<%s(start=%d end=%d strand=%s node_type=%d>" %
                (self.__class__.__name__, self.start, self.end, self.strand, self.node_type))
    def __eq__(self, other):
        return ((self.start == other.start) and
                (self.end == other.end) and
                (self.strand == other.strand) and
                (self.node_type == other.node_type))
    def __ne__(self, other):
        return not self.__eq__(other)
    def __hash__(self):
        return (self.start << 18) | (self.end << 2) | (self.strand << 1) | (self.node_type)

# isoform graph node
NodeTuple = collections.namedtuple('NodeTuple', ['start', 'end', 'strand', 'node_type'])
# information for merging transcripts
MergeTuple = collections.namedtuple('MergeTuple', ['node', 'scores', 'predecessors', 'successors'])

def interval_overlap(a, b):
    return (a[0] < b[1]) and (b[0] < a[1])

def cmp_strand(a, b):
    if a == b:
        return True
    if (a == ".") or (b == "."):
        return True
    return False

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
    nodes = leaf_nodes + internal_nodes    
    # compare each leaf to rest of nodes
    for i in xrange(len(leaf_nodes)):
        n1 = nodes[i]                
        n1_introns = G.predecessors(n1) + G.successors(n1)        
        for j in xrange(i+1, len(nodes)):
            n2 = nodes[j]
            # ensure proxy nodes are not considered for merging
            if (n2 in node_proxies[n1]) or (n1 in node_proxies[n2]):
                continue
            n2_introns = G.predecessors(n2) + G.successors(n2)
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
    elif strand1 == ".":
        return strand2
    elif strand2 == ".":
        return strand1
    else:
        logging.error("Incompatible strands")        
        assert False

def merge_intron_compatible_exons(G, nodes):    
    logging.debug("Nodes to merge: %s" % nodes)    
    newnode = NodeTuple(min([n.start for n in nodes]),
                        max([n.end for n in nodes]),
                        reduce(merge_strand, [n.strand for n in nodes]),
                        EXON)
    newscores = []
    predecessors = set()
    successors = set()
    for n in nodes:
        newscores.extend(G.node[n]['scores'])
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
    if node.strand == "-":
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
        newnode = NodeTuple(start, end, merge_strand(leaf_node.strand, match_node.strand), EXON)
        return newnode    
    else:
        return None 

def find_trimmable_nodes(G, leaf_nodes, internal_nodes, overhang_threshold):
    # compare each leaf to rest of nodes
    nodes = leaf_nodes + internal_nodes    
    for i in xrange(len(leaf_nodes)):
        n1 = nodes[i]        
        # compare each node against all other nodes in case
        # it could trim against more than one node at a time.
        # need to find all possible trimmings and choose a 
        # consensus otherwise strange exponential behavior occurs
        # in special cases during collapse
        trimmed_nodes = []
        logging.debug("leaf=%s" % str(n1))
        for j in xrange(len(nodes)):
            # do not compare against self
            if i == j:
                continue
            n2 = nodes[j]           
            logging.debug("match=%s" % str(n2))
            # see if leaf node can be trimmed to become 
            # intron compatible with another node
            newnode = trim(G, n1, n2, overhang_threshold)
            if newnode is not None:
                trimmed_nodes.append(newnode)
                logging.debug("trimmed node=%s" % str(newnode))
        
        if len(trimmed_nodes) > 0:
            # build a consensus trimmed node from all possible trimmings
            # by choosing coordinates that give the smallest new exon
            start = max(n.start for n in trimmed_nodes)
            end = min(n.end for n in trimmed_nodes)
            scores = G.node[n1]['scores']        
            newnode = NodeTuple(start, end, trimmed_nodes[0].strand, EXON)
            # check if the new node is already part of the graph
            if newnode in G:
                logging.debug("Trimmed node already part of graph=%s" % str(newnode))
                continue
            logging.debug("Final trimmed node=%s" % str(newnode))
            merge_tuple = MergeTuple(newnode, scores, G.predecessors(n1), G.successors(n1))
            yield merge_tuple, n1


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
            left_free = (((exon_index == 0) and (strand != "-")) or
                         ((exon_index == len(exons)) and (strand == "-")))
            right_free = (((exon_index == len(exons)) and (strand != "-")) or
                          ((exon_index == 0) and (strand == "-")))
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
        partial_path_dict[p[0]].append(p)
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
        scores = G.node[node]['scores']
        total_score += sum(scores)/len(scores)
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
        if n in self.G:
            self.G.node[n]['scores'].extend(list(scores))
        else:
            self.G.add_node(n, scores=list(scores))
        return n

    def add_transcript(self, transcript):
        exons = transcript.exons
        strand = transcript.strand
        # must have at least one exon
        e1_node = NodeTuple(exons[0].start, exons[0].end, strand, EXON)
        self._add_node(e1_node, [exons[0].score])
        for e2 in exons[1:]:
            # add intron node            
            intron = NodeTuple(e1_node.end, e2.start, strand, INTRON)
            self._add_node(intron, [transcript.score]) 
            # add exon
            e2_node = NodeTuple(e2.start, e2.end, strand, EXON) 
            self._add_node(e2_node, [e2.score])
            # add edges from exon -> intron -> exon
            if transcript.strand == "-":
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
    
    def _collapse_once(self, trim=False, overhang_threshold=0):
        # build interval tree
        exon_tree = IntervalTree()        
        for node in self.G.nodes_iter():
            if node.node_type == EXON:
                exon_tree.insert_interval(Interval(node.start, node.end, value=node))

        leaf_node_indexes = set()
        nodes = []
        # cluster the exons by strand
        cluster_trees = {"+": ClusterTree(0,1), "-": ClusterTree(0,1)}
        exon_index = 0
        for node in self.G.nodes_iter():
            if node.node_type == INTRON:
                continue
            if self.G.in_degree(node) == 0 or self.G.out_degree(node) == 0:
                leaf_node_indexes.add(exon_index)
            # add to cluster tree based on strand
            if node.strand == ".":
                strands = ["+", "-"]
            else:
                strands = [node.strand]
            for strand in strands:
                cluster_trees[strand].insert(node.start, node.end, exon_index)
            nodes.append(node)
            exon_index += 1
        
        logging.debug('COLLAPSE')
        merge_tuples = []
        replaced_nodes = set()                
        for strand, cluster_tree in cluster_trees.iteritems():
            logging.debug("strand=%s" % strand)
            # find clusters of overlapping nodes
            for start, end, indexes in cluster_tree.getregions():
                leaf_nodes = [nodes[i] for i in set(indexes).intersection(leaf_node_indexes)]                
                logging.debug("leaf_nodes=%s" % leaf_nodes)
                if len(leaf_nodes) == 0:
                    continue
                internal_nodes = [nodes[i] for i in set(indexes).difference(leaf_node_indexes)]
                logging.debug("internal_nodes=%s" % internal_nodes)
                if len(leaf_nodes) + len(internal_nodes) < 2:
                    continue

                node_proxies = collections.defaultdict(lambda: set())
                # only perform trimming if enabled
                if trim:
                    # find overhanging exons that can be trimmed from graph
                    logging.debug("TRIMMING")
                    for trim_merge_tuple, leaf_node in find_trimmable_nodes(self.G, leaf_nodes, internal_nodes, overhang_threshold):
                        # add the trimmable nodes to the graph temporarily so 
                        # that they can participate in merging
                        self._add_merged_node(trim_merge_tuple)
                        # all the trimmed nodes will be removed unless they are 
                        # used to facilitate a merge
                        replaced_nodes.add(trim_merge_tuple.node)
                        # the trimmable nodes are by definition leaf nodes
                        leaf_nodes.append(trim_merge_tuple.node)
                        # maintain a link from the trimmed node to the original 
                        # leaf nodes so that all can be removed if necessary
                        node_proxies[trim_merge_tuple.node].add(leaf_node)
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
                        replaced_nodes.update(node_proxies[n])
                logging.debug("/MERGING")
        logging.debug("/COLLAPSE")
        # remove all exons that were replaced by merges
        self.G.remove_nodes_from(replaced_nodes)
        # reconstruct the graph by adding all the merged exons
        for merge_tuple in merge_tuples:
            self._add_merged_node(merge_tuple)
        return len(merge_tuples)

    def collapse(self, trim=False, overhang_threshold=0):        
        merges = 1
        while merges > 0:
            merges = self._collapse_once(trim=trim, overhang_threshold=overhang_threshold)
            logging.debug("merges performed=%d" % merges)
#    def collapse(self, trim=False, overhang_threshold=0):        
#        trims = 1
#        while trims > 0:
#            # merge graph so only intron-incompatible transcripts remain
#            self._merge(trim=trim, overhang_threshold=overhang_threshold)        
#            if trim:
#                # trim overhanging exons from graph
#                trims = self._collapse_once(mode="trim", overhang_threshold=overhang_threshold)
#                logging.debug("trims performed=%d" % trims)
#            else:
#                trims = 0



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
            # flip exons for minus strand to keep path in correct order
            exons = t.exons[::-1] if t.strand == "-" else t.exons
            for path in thread_transcript(self.G, exon_tree, t.strand, exons):
                # ignore paths of length 1 (useless for threading)
                if len(path) > 1:
                    # TODO: check path for now
                    for i in xrange(len(path)-1):
                        for j in xrange(i+1, len(path)):
                            if nx.shortest_path(self.G, path[i], path[j]) == False:
                                logging.error("PATH ERROR")
                                assert False
                    threaded_paths.add(tuple(path))

        logging.debug("threaded %s paths" % len(threaded_paths))
        for p in threaded_paths:
            logging.debug("threaded path=%s" % str([str(n) for n in p]))
            logging.debug("threaded exon sizes=%s" % str([str(n.end - n.start) for n in p]))
            logging.debug("threaded exon starts=%s" % str([str(n.start) for n in p]))
            
        logging.debug("/THREADING")        
        # using transcript-threaded paths as a guide find full paths
        # through the graph        
        #for paths in find_paths(self.G, []):
        all_paths = set()
        for paths in find_paths(self.G, threaded_paths):
            for p in paths:
                all_paths.add(tuple(p))
        for p in all_paths:
            yield p, score_path(self.G, p)

#    class Node(object):
#        __slots__ = ('start', 'end', 'strand', 'node_type', 'scores')
#        def __init__(self, start, end, strand, node_type):
#            self.start = start
#            self.end = end
#            self.strand = strand
#            self.node_type = node_type
#            self.scores = []
#        def __eq__(self, other):
#            return ((self.start == other.start) and
#                    (self.end == other.end) and
#                    (self.strand == other.strand) and
#                    (self.node_type == other.node_type))
#        def __ne__(self, other):
#            return not self.__eq__(other)
#        def __hash__(self):
#            return (self.start << 18) | (self.end << 2) | (self.strand << 1) | (self.node_type)
