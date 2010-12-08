'''
Created on Nov 23, 2010

@author: mkiyer
'''
import logging
import networkx as nx
import collections
import itertools
import operator

from bx.intervals.intersection import Interval, IntervalTree
from bx.intervals.cluster import ClusterTree

from cNode import Node, strand_int_to_str, strand_str_to_int, node_type_to_str
from base import EXON, INTRON, POS_STRAND, NEG_STRAND, NO_STRAND

# a structure for managing isoform graph information
# stores the node plus its immediate edges
ExonTuple = collections.namedtuple('ExonTuple', ['node', 'preds', 'succs', 'exon_domains', 'intron_domains'])
MergeTuple = collections.namedtuple('MergeTuple', ['node', 'ids', 'preds', 'succs'])

class ClusterInfo(object):
    __slots__ = ('indexes', 'preds', 'succs')
    def __init__(self, indexes=None, preds=None, succs=None):
        self.indexes = set()
        self.preds = set()
        self.succs = set()
        if indexes is not None:
            self.indexes.update(indexes)
        if preds is not None:
            self.preds.update(preds)
        if succs is not None:
            self.succs.update(succs)
    def update(self, other):
        self.indexes.update(other.indexes)
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

def merge_exon_tuples(G, cluster_start, cluster_end, strand, exon_tuples, join_preds=None, join_succs=None):
    preds = set()
    if join_preds is not None:
        preds.update(join_preds)
    succs = set()
    if join_succs is not None:
        succs.update(join_succs)
    ids = set()
    # TODO: overkill but ensures all strands are correct
    strand = reduce(merge_strand, set([NO_STRAND]).union([x.strand for x in preds], [x.strand for x in succs]))
    #strand = NO_STRAND
    for exon_tuple in exon_tuples:
        preds.update(exon_tuple.preds)
        succs.update(exon_tuple.succs)
        node = exon_tuple.node
        ids.update(G.node[node]['ids'])
        strand = merge_strand(strand, node.strand)
    return MergeTuple(Node(cluster_start, cluster_end, strand, EXON), ids, preds, succs)    

def get_free_ends(strand, in_degree, out_degree):
    """returns whether the (left,right) genomic ends of the node are free ends without edges"""
    # TODO: move to C code
    fiveprime = (in_degree == 0)
    threeprime = (out_degree == 0)
    if strand == NEG_STRAND:
        return threeprime, fiveprime
    else:
        return fiveprime, threeprime

def trim_interval(start, end, trim_left, trim_right, overhang_threshold):
    # TODO: move to C code
    trim_start, trim_end = start, end
    if trim_left:
        trim_start = start + overhang_threshold
        if trim_start >= end:
            trim_start = end-1
    if trim_right:
        trim_end = end - overhang_threshold
        if trim_end <= start:
            trim_end = start+1
    if trim_end < trim_start:
        return trim_end, trim_start
    return trim_start, trim_end

def partition_nodes_by_strand(G):
    # TODO: move to C code
    strand_introns = [[],[]]
    strand_exons = [[], [], []]
    for node in G:
        if node.node_type == INTRON:
            assert node.strand != NO_STRAND
            strand_introns[node.strand].append(node)
        else:
            strand_exons[node.strand].append(node)
    return strand_exons, strand_introns

def make_intron_boundary_map(introns):
    # build an index of intron edges so that 
    # free exon ends can be joined with compatible introns
    intron_starts = collections.defaultdict(lambda: set())
    intron_ends = collections.defaultdict(lambda: set())
    for node in introns:
        # no unstranded introns allowed
        assert node.strand != NO_STRAND
        # keep track of positions where introns can be joined to exons
        intron_starts[node.start].add(node)
        intron_ends[node.end].add(node)
    return dict(intron_starts), dict(intron_ends)

def find_intron_domains(strand, cluster_start, cluster_end, introns):
    '''
    define genomic intervals representing unique intron 'domains'
    optimize the merging by first clustering transcripts within intron compatible "domains"
    '''
    # add all intron start/end positions to dictionaries
    intron_boundaries = set()
    for intron in introns:
        # check strand
        assert intron.strand != NO_STRAND
        assert intron.strand == strand
        # keep track of positions where introns can be joined to exons
        intron_boundaries.add(intron.start)
        intron_boundaries.add(intron.end)
    # sort the intron boundary positions and add them to interval trees
    sorted_intron_boundaries = [cluster_start]
    sorted_intron_boundaries.extend(sorted(intron_boundaries))
    sorted_intron_boundaries.append(cluster_end)
    tree = IntervalTree()
    intervals = []
    domain_id = 0
    for start, end in itertools.izip(sorted_intron_boundaries[:-1], sorted_intron_boundaries[1:]):
        intervals.append(Interval(start, end, strand=strand, value=domain_id))        
        tree.insert_interval(intervals[-1])
        domain_id += 1
    return intervals, tree

def get_node_intron_domains(exon, preds, succs, intron_domain_tree, overhang_threshold):
    # trim exons before computing exon domains
    if overhang_threshold > 0:
        left_free, right_free = get_free_ends(exon.strand, len(preds), len(succs))
        start, end = trim_interval(exon.start, exon.end, 
                                   left_free, right_free,
                                   overhang_threshold)
    else:
        start, end = exon.start, exon.end
    # transform exon coordinates to intron domain regions
    exon_domains = set(interval.value for interval in intron_domain_tree.find(start, end))
    logging.debug('exon start=%d end=%d trimstart=%d trimend=%d domains=%s' % (exon.start, exon.end, start, end, exon_domains))
    # convert introns into domain regions
    intron_domains = set()
    for intron in itertools.chain(preds, succs):
        intron_domains.update(interval.value for interval in intron_domain_tree.find(intron.start, intron.end))
    return exon_domains, intron_domains

def map_exons_to_intron_domains(G, exons, domain_tree, overhang_threshold=0):
    for node in exons:
        assert node.node_type == EXON
        preds, succs = G.predecessors(node), G.successors(node)
        edomains, idomains = get_node_intron_domains(node, preds, succs,
                                                     domain_tree, 
                                                     overhang_threshold)
        yield ExonTuple(node, preds, succs, edomains, idomains)

def cluster_within_domains(exon_tuples, domain_set):
    # TODO: can probably speed up this check by pre-storing
    # all the single domain -> mappings and then taking the union
    # for ranges of domains
    matched_indexes = []
    for i,exon_tuple in enumerate(exon_tuples):
        if (domain_set.issuperset(exon_tuple.exon_domains) and
            domain_set.isdisjoint(exon_tuple.intron_domains)):
            matched_indexes.append(i)
    # the exons within this domain range can be
    # clustered by their genomic positions 
    cluster_tree = ClusterTree(0,1)
    for i in matched_indexes:
        node = exon_tuples[i].node
        cluster_tree.insert(node.start, node.end, i)
    # find groups of overlapping nodes
    for start, end, overlapping_indexes in cluster_tree.getregions():
        yield start, end, overlapping_indexes 

#def cluster_exons(G, strand, exons, introns, overhang_threshold):
#    cluster_start = min(e.start for e in exons)
#    cluster_end = max(e.end for e in exons)
#    # TODO: strand should not needed here?
#    domain_intervals, domain_trees = find_intron_domains(strand, cluster_start, cluster_end, introns)
#    intron_start_dict, intron_end_dict = make_intron_boundary_map(introns)  
#    # get 'intron compatibility' domains associated with each exon
#    exon_tuples = list(map_exons_to_intron_domains(G, exons, domain_trees, overhang_threshold))    
#    # find set of possible domain-ranges to be clustered
#    #domain_sets = set(frozenset(e.exon_domains) for e in exon_tuples)
#    # sorted sets by their size (smallest to largest)
#    #sorted_domain_sets = sorted(domain_sets, key=len)
#    #for domain_set in sorted_domain_sets:
#    #    print 'domain set %s' % str(domain_set)
#    #    cluster_domains(exon_tuples, domain_set, clusters)
#    clusters = collections.defaultdict(lambda: ClusterInfo())
#    # start with single domain regions and grow
#    for domain_size in xrange(0, len(domain_intervals)):
#        # scan through all domains of this size and cluster exons
#        for domain_end in xrange(domain_size, len(domain_intervals)):
#            domain_start = domain_end - domain_size
#            domain_set = set(range(domain_start, domain_end+1))
#            # find genomic intervals of domain            
#            interval_start = domain_intervals[domain_start].start          
#            interval_end = domain_intervals[domain_end].end
#            interval_left_introns = intron_end_dict.get(interval_start, set()) 
#            interval_right_introns = intron_start_dict.get(interval_end, set())
#            logging.debug('domain set (%d,%d) interval (%d,%d)' % (domain_start, domain_end, interval_start, interval_end))            
#            # cluster exons within this domain
#            for start, end, indexes in cluster_within_domains(exon_tuples, domain_set):
#                # trim clusters that lie outside of domain interval
#                if start < interval_start:
#                    #print 'trim start %d < interval start %d' % (start, interval_start)
#                    start = interval_start
#                connect_left = (start == interval_start) and len(interval_left_introns) > 0
#                if end > interval_end:
#                    #print 'trim end %d > interval end %d' % (end, interval_end)
#                    end = interval_end
#                connect_right = (end == interval_end) and len(interval_right_introns) > 0
#                if (len(indexes) > 1 or connect_left or connect_right):
#                    # TODO: delete print statements
#                    #print 'cluster start=%d end=%d indexes=%s' % (start, end, indexes)
#                    #for i in indexes:
#                    #    print '  node=%s' % (str(exon_tuples[i]))
#                    # check whether cluster exactly touches the domain borders
#                    # to determine whether to join the cluster with adjacent
#                    # introns
#                    left_introns, right_introns = set(), set()
#                    if start == interval_start:
#                        left_introns.update(interval_left_introns)
#                    if end == interval_end:
#                        right_introns.update(interval_right_introns)
#                    if strand == NEG_STRAND:
#                        preds, succs = right_introns, left_introns
#                    else:
#                        preds, succs = left_introns, right_introns
#                    # add to cluster dictionary
#                    cluster_info = ClusterInfo(indexes, preds, succs)
#                    logging.debug('cluster start=%d end=%d indexes=%s preds=%s succs=%s' % (start, end, indexes, preds, succs))                
#                    clusters[(start,end)].update(cluster_info) 
#    # build merge tuples from the clusters
#    merge_list = []
#    for startend, cluster_info in clusters.iteritems():
#        start, end = startend
#        orig_exon_tuples = [exon_tuples[i] for i in cluster_info.indexes]
#        orig_nodes = set(e.node for e in orig_exon_tuples)        
#        merge_tuple = merge_exon_tuples(G, start, end, strand, orig_exon_tuples, 
#                                        cluster_info.preds, cluster_info.succs) 
#        logging.debug('start=%d end=%d merged node=%s' % (start, end, str(merge_tuple)))
#        # decorate merge tuples with number of nodes merged to
#        # sort later
#        nodeset = orig_nodes.union(merge_tuple.preds, merge_tuple.succs)
#        merge_list.append((len(nodeset), nodeset, merge_tuple, orig_nodes))
#    # sort by number of nodes merged from largest to smallest
#    merge_list = [x[1:] for x in sorted(merge_list, key=operator.itemgetter(0), reverse=True)]
#    indexes_to_discard = set()
#    # remove subsets from list of merge tuples
#    for i in xrange(len(merge_list) - 1):
#        nodeset1 = merge_list[i][0]
#        for j in xrange(i+1, len(merge_list)):
#            nodeset2 = merge_list[j][0]
#            if nodeset1.issuperset(nodeset2):
#                logging.debug('discarding subset tuple %s' % str(merge_list[j][1]))
#                indexes_to_discard.add(j)
#    merge_list = [merge_list[i][1:] for i in xrange(len(merge_list))
#                  if i not in indexes_to_discard] 
#    return merge_list

def collapse_cluster(G, strand, exon_tuples, 
                     domain_indexes, domain_intervals,
                     intron_start_dict, intron_end_dict):
    clusters = collections.defaultdict(lambda: ClusterInfo())
    # start with single domain regions and grow
    for domain_size in xrange(0, len(domain_indexes)):
        # scan through all domains of this size and cluster exons
        for end_index in xrange(domain_size, len(domain_indexes)):
            start_index = end_index - domain_size
            domain_start = domain_indexes[start_index]
            domain_end = domain_indexes[end_index]
            domain_set = set(range(domain_start, domain_end+1))
            # find genomic intervals of domain            
            interval_start = domain_intervals[domain_start].start          
            interval_end = domain_intervals[domain_end].end
            interval_left_introns = intron_end_dict.get(interval_start, set()) 
            interval_right_introns = intron_start_dict.get(interval_end, set())
            logging.debug('domain set (%d,%d) interval (%d,%d)' % (domain_start, domain_end, interval_start, interval_end))            
            # cluster exons within this domain
            for start, end, indexes in cluster_within_domains(exon_tuples, domain_set):
                # trim clusters that lie outside of domain interval
                if start < interval_start:
                    #print 'trim start %d < interval start %d' % (start, interval_start)
                    start = interval_start
                connect_left = (start == interval_start) and len(interval_left_introns) > 0
                if end > interval_end:
                    #print 'trim end %d > interval end %d' % (end, interval_end)
                    end = interval_end
                connect_right = (end == interval_end) and len(interval_right_introns) > 0
                if (len(indexes) > 1 or connect_left or connect_right):
                    # TODO: delete print statements
                    #print 'cluster start=%d end=%d indexes=%s' % (start, end, indexes)
                    #for i in indexes:
                    #    print '  node=%s' % (str(exon_tuples[i]))
                    # check whether cluster exactly touches the domain borders
                    # to determine whether to join the cluster with adjacent
                    # introns
                    left_introns, right_introns = set(), set()
                    if start == interval_start:
                        left_introns.update(interval_left_introns)
                    if end == interval_end:
                        right_introns.update(interval_right_introns)
                    if strand == NEG_STRAND:
                        preds, succs = right_introns, left_introns
                    else:
                        preds, succs = left_introns, right_introns
                    # add to cluster dictionary
                    cluster_info = ClusterInfo(indexes, preds, succs)
                    logging.debug('cluster start=%d end=%d indexes=%s preds=%s succs=%s' % (start, end, indexes, preds, succs))                
                    clusters[(start,end)].update(cluster_info) 
    # build merge tuples from the clusters
    merge_list = []
    for startend, cluster_info in clusters.iteritems():
        start, end = startend
        orig_exon_tuples = [exon_tuples[i] for i in cluster_info.indexes]
        orig_nodes = set(e.node for e in orig_exon_tuples)        
        merge_tuple = merge_exon_tuples(G, start, end, strand, orig_exon_tuples, 
                                        cluster_info.preds, cluster_info.succs) 
        logging.debug('start=%d end=%d merged node=%s' % (start, end, str(merge_tuple)))
        # decorate merge tuples with number of nodes merged to
        # sort later
        nodeset = orig_nodes.union(merge_tuple.preds, merge_tuple.succs)
        merge_list.append((len(nodeset), nodeset, merge_tuple, orig_nodes))
    # sort by number of nodes merged from largest to smallest
    merge_list = [x[1:] for x in sorted(merge_list, key=operator.itemgetter(0), reverse=True)]
    indexes_to_discard = set()
    # remove subsets from list of merge tuples
    for i in xrange(len(merge_list) - 1):
        nodeset1 = merge_list[i][0]
        for j in xrange(i+1, len(merge_list)):
            nodeset2 = merge_list[j][0]
            if nodeset1.issuperset(nodeset2):
                logging.debug('discarding subset tuple %s' % str(merge_list[j][1]))
                indexes_to_discard.add(j)
    merge_list = [merge_list[i][1:] for i in xrange(len(merge_list))
                  if i not in indexes_to_discard]
    return merge_list

def collapse_graph(G, strand, exons, introns, overhang_threshold):
    # divide graph into regions of intron compatibility,
    # or "intron domains"
    cluster_start = min(e.start for e in exons)
    cluster_end = max(e.end for e in exons)
    domain_intervals, domain_tree = find_intron_domains(strand, cluster_start, cluster_end, introns)
    intron_start_dict, intron_end_dict = make_intron_boundary_map(introns)
    # partition the exons into overlapping clusters
    # so that each cluster can be solved independently
    exon_cluster_tree = ClusterTree(0,1)
    for exon_index,exon in enumerate(exons):
        exon_cluster_tree.insert(exon.start, exon.end, exon_index)
    # process clusters of overlapping exon nodes
    merge_list = []
    for start, end, indexes in exon_cluster_tree.getregions():
        # get 'intron compatibility' domains associated with this entire exon cluster
        cluster_domains = sorted(interval.value for interval in domain_tree.find(start, end))
        # get 'intron compatibility' domains associated with individual exons in cluster
        exon_tuples = list(map_exons_to_intron_domains(G, [exons[i] for i in indexes], domain_tree, overhang_threshold))
        merge_list.extend(collapse_cluster(G, strand, exon_tuples, 
                                           cluster_domains, domain_intervals,
                                           intron_start_dict, intron_end_dict))
    return merge_list

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

    def collapse(self, trim=False, overhang_threshold=0):
        logging.info("COLLAPSE called with nodes=%d" % (len(self.G)))
        # TODO: just use trim and not overhang threshold
        if not trim:
            overhang_threshold = 0
        # divide nodes into pos/neg strands
        strand_exons, strand_introns = partition_nodes_by_strand(self.G)        
        merge_tuples = []
        replaced_nodes = set()         
        for strand in (POS_STRAND, NEG_STRAND):
            exons = strand_exons[strand]
            introns = strand_introns[strand]
            # by default cluster unstranded exons as though
            # they are on the (+) strand, so if there are
            # no (-) strand exons can skip clustering
            if len(exons) == 0 and strand == NEG_STRAND:
                continue
            # add unstranded exons
            exons.extend(strand_exons[NO_STRAND])        
            if len(exons) == 0:
                continue
            # run the main clustering algorithm
            for merge_tuple, orig_nodes in collapse_graph(self.G, strand, exons, introns, overhang_threshold):
                # if merged node is already in the graph, then simply
                # update it with the new list of ids rather than replacing
                # it and adding it again later
                if merge_tuple.node in self.G:
                    self._add_merged_node(merge_tuple)
                    replaced_nodes.update(set(orig_nodes).difference([merge_tuple.node]))
                else:
                    merge_tuples.append(merge_tuple)
                    replaced_nodes.update(orig_nodes)
#            for merge_tuple, orig_nodes in cluster_exons(self.G, strand, exons, introns, overhang_threshold):
#                # if merged node is already in the graph, then simply
#                # update it with the new list of ids rather than replacing
#                # it and adding it again later
#                if merge_tuple.node in self.G:
#                    self._add_merged_node(merge_tuple)
#                    replaced_nodes.update(set(orig_nodes).difference([merge_tuple.node]))
#                else:
#                    merge_tuples.append(merge_tuple)
#                    replaced_nodes.update(orig_nodes)
        self.G.remove_nodes_from(replaced_nodes)
        for merge_tuple in merge_tuples:
            self._add_merged_node(merge_tuple)
        logging.debug("/COLLAPSE")
