'''
Created on May 5, 2010

@author: mkiyer
'''
import collections
import logging
import operator
import sys

if sys.version_info < (2,6):
    from namedtuple import namedtuple
else:
    from collections import namedtuple

import networkx as nx
import matplotlib.pyplot as plt
import numpy as np

from bx.intervals.intersection import Interval, IntervalTree
from bx.intervals.cluster import ClusterTree

# local imports
import cuffcompare
from transcript import Transcript
from transcript_parser import parse_gtf_chrom_loci

# node tuple
NodeTuple = namedtuple('NodeTuple', ['start', 'end', 'strand', 'is_exon'])
# information for merging transcripts
MergeTuple = namedtuple('MergeTuple', ['node', 'score', 'predecessors', 'successors'])

# graph score attribute name (node dict)
GRAPH_SCORE_KEY = 'weight'

def get_start_node():
    return NodeTuple(start=-1, end=-1, strand=NO_STRAND, is_exon=False)
def get_end_node():
    return NodeTuple(start=-2, end=-2, strand=NO_STRAND, is_exon=False)

# convert strand to integer
POS_STRAND = 1
NEG_STRAND = -1
NO_STRAND = 0
_encode_strand = {'+': POS_STRAND,
                  '-': NEG_STRAND,
                  '.': NO_STRAND}
_decode_strand = {POS_STRAND: '+',
                  NEG_STRAND: '-',
                  NO_STRAND: '.'}

def encode_strand(strand_string):
    global _encode_strand
    if strand_string == None:
        return 0
    else:
        return _encode_strand[strand_string]

def decode_strand(strand):
    global _decode_strand
    return _decode_strand[strand]

def cmp_strand(a, b):
    if a == b:
        return True
    if (a == NO_STRAND) or (b == NO_STRAND):
        return True
    return False

def interval_overlap(a, b):
    return (a[0] < b[1]) and (b[0] < a[1])

def cluster_single_exon_transcripts(transcripts):
    '''
    optimize the merging by first clustering single exon transcripts
    that do not overlap any multi-exon transcripts
    '''
    logging.debug("Before clustering transcripts=%d" % len(transcripts))
    newtranscripts = []
    single_exon_transcripts = []
    # build intron interval tree
    intron_tree = IntervalTree()
    for t in transcripts:
        for intron in t.introns:
            intron_tree.insert_interval(Interval(intron[0], intron[1], strand=t.strand))            
        if (len(t.exons) > 1) or (t.strand != '.'):
            # do not manipulate multi-exon transcripts here
            newtranscripts.append(t)
        else:
            single_exon_transcripts.append(t)
    
    # organize exons by how they overlap introns by 
    # building an (intron list) -> (single exon transcript) mapping
    intronlist_exon_map = collections.defaultdict(lambda: [])
    for t in single_exon_transcripts:
        overlapping_introns = intron_tree.find(t.tx_start, t.tx_end)
        intron_tuple = tuple((intron.start, intron.end, intron.strand)
                             for intron in overlapping_introns)
        intronlist_exon_map[intron_tuple].append(t)
    
    # transcripts that overlap identical intron lists
    # can be directly clustered without sacrificing 
    # potential merging with multi-exon transcripts
    chrom = transcripts[0].chrom
    for intronlist, tlist in intronlist_exon_map.iteritems():
        #logging.debug("intron list: %s has %d transcripts" % (str(intronlist), len(tlist)))
        #logging.debug("exons to cluster: %s" % (['(%s,%s)' % (x.tx_start, x.tx_end) for x in tlist]))
        cluster_tree = ClusterTree(0, 1)
        for index, t in enumerate(tlist):
            cluster_tree.insert(t.tx_start, t.tx_end, index)
        for start, end, tlist_indexes in cluster_tree.getregions():
            # create a new Transcript
            old_transcripts = [tlist[i] for i in tlist_indexes]
            # average score
            score = max(t.score for t in old_transcripts)
            #score = sum(t.score for t in old_transcripts) / float(len(old_transcripts))            
            t = Transcript(id=old_transcripts[0].id,
                           chrom=chrom,
                           strand='.',
                           score=score,
                           exons=[(start, end)])
            newtranscripts.append(t)
    logging.debug("After clustering transcripts=%d" % (len(newtranscripts)))
    return newtranscripts


def merge_score(score1, score2):
    return max(score1, score2)

def add_node_update_score(G, n, score):
    G.add_node(n)
    oldscore = G.node[n].setdefault(GRAPH_SCORE_KEY, score)
    #G.node[n][GRAPH_SCORE_KEY] = update_node_score(n, oldweight, n, score)
    G.node[n][GRAPH_SCORE_KEY] = merge_score(oldscore, score)

def build_isoform_graph(transcripts):
    G = nx.DiGraph()
    for t in transcripts:
        # get strand
        strand = encode_strand(t.strand)
        # nodes are tuples (start, end, strand, is_exon)
        e1 = NodeTuple(t.exons[0][0], t.exons[0][1], strand, True)
        add_node_update_score(G, e1, t.score)
        for e2 in t.exons[1:]:
            # add exon/intron nodes
            intron = NodeTuple(e1[1], e2[0], strand, False)
            # TODO: get junction expression levels!
            add_node_update_score(G, intron, 0.0)
            e2 = NodeTuple(e2[0], e2[1], strand, True)
            add_node_update_score(G, e2, t.score)
            # add edges from exon -> intron -> exon
            G.add_edge(e1, intron)
            G.add_edge(intron, e2)
            # continue loop
            e1 = e2
    return G

def is_intron_compatible(exon1, introns1, exon2, introns2):
    for intron in introns2:
        if interval_overlap(exon1, intron):
            return False
    for intron in introns1:
        if interval_overlap(exon2, intron):
            return False
    return True

def get_overhang(leaf_exon, match_exon):
    # check if exon is overhanging by a small amount
    # and trim if necessary
    start_diff = match_exon.start - leaf_exon.start
    end_diff = leaf_exon.end - match_exon.end       
    if match_exon.strand == NEG_STRAND:
        return end_diff, start_diff
    else:
        return start_diff, end_diff

def make_merge_tuple(G, node1, score1, node2, score2):
    # make new exon
    newnode = NodeTuple(min(node1[0], node2[0]),
                        max(node1[1], node2[1]),
                        node1[2] or node2[2],
                        True)
    newscore = merge_score(score1, score2)
    #newscore = update_node_score(node1, score1, node2, score2)
    # connect all neighbors to this new exon
    predecessors = G.predecessors(node1) + G.predecessors(node2)
    successors = G.successors(node1) + G.successors(node2)
    return MergeTuple(newnode, newscore, predecessors, successors)

def add_merged_exon(G, merge_data):
    newnode, newscore, predecessors, successors = merge_data
    # add node to graph
    add_node_update_score(G, newnode, newscore)
    # connect all neighbors to this new exon
    for neighbor in predecessors:
        G.add_edge(neighbor, newnode)
    for neighbor in successors:
        G.add_edge(newnode, neighbor)

def merge_intron_compatible_exons(G, leaf_exon, match_exon):
    leaf_introns = G.predecessors(leaf_exon) + G.successors(leaf_exon)
    leaf_score = G.node[leaf_exon][GRAPH_SCORE_KEY] 
    match_introns = G.predecessors(match_exon) + G.successors(match_exon)
    match_score = G.node[match_exon][GRAPH_SCORE_KEY]     
    # leaf exon cannot overlap match introns, and 
    # match exon cannot overlap leaf introns
    if is_intron_compatible(leaf_exon, leaf_introns, match_exon, match_introns):
        merge_tuple = make_merge_tuple(G, leaf_exon, leaf_score, match_exon, match_score)    
        return merge_tuple, [leaf_exon, match_exon]
    return None, None

def make_trimmed_merge_tuple(G, exon1, score1, exon2, score2, 
                             fiveprime=True, threeprime=True):
    new_start, new_end = exon1.start, exon1.end
    predecessors = []
    successors = []
    if exon2.strand == NEG_STRAND:
        if fiveprime:
            new_end = exon2.end
            successors = G.successors(exon1) + G.successors(exon2)
        if threeprime:
            new_start = exon2.start      
            predecessors = G.predecessors(exon1) + G.predecessors(exon2)
    else:
        if fiveprime:
            new_start = exon2.start
            predecessors = G.predecessors(exon1) + G.predecessors(exon2)
        if threeprime:
            new_end = exon2.end
            successors = G.successors(exon1) + G.successors(exon2)
    # make new exon
    newnode = NodeTuple(new_start, new_end, exon1.strand or exon2.strand, True)
    newscore = merge_score(score1, score2)
    # connect all neighbors to this new exon
    return MergeTuple(newnode, newscore, predecessors, successors)


def merge_overhanging_exons(G, leaf_exon, match_exon, overhang_threshold=0):
    leaf_score = G.node[leaf_exon][GRAPH_SCORE_KEY] 
    match_score = G.node[match_exon][GRAPH_SCORE_KEY]
    introns_5prime = G.predecessors(leaf_exon)
    introns_3prime = G.successors(leaf_exon)
    overhang_5prime, overhang_3prime = get_overhang(leaf_exon, match_exon)
    canmerge_5prime = (overhang_5prime >= 0) and (overhang_5prime < overhang_threshold)
    canmerge_3prime = (overhang_3prime >= 0) and (overhang_3prime < overhang_threshold)

    # ensure transcripts are not intron compatible
    # and if they are perform normal merging    
    leaf_introns = G.predecessors(leaf_exon) + G.successors(leaf_exon)
    match_introns = G.predecessors(match_exon) + G.successors(match_exon)    
    if is_intron_compatible(leaf_exon, leaf_introns, match_exon, match_introns):
        return None, None

    merge_5prime, merge_3prime = False, False
    if len(introns_5prime) + len(introns_3prime) == 0:
        # single exon "transcript", check both ends
        merge_5prime = canmerge_5prime
        merge_3prime = canmerge_3prime        
    elif len(introns_5prime) == 0:
        # this is a putative TSS (5' start site)
        # introns must be a subset of match introns 
        match_introns_3prime = G.successors(match_exon)
        merge_5prime = ((set(introns_3prime) == set(match_introns_3prime)) and canmerge_5prime)
    else:
        match_introns_5prime = G.predecessors(match_exon)
        merge_3prime = ((set(introns_5prime) == set(match_introns_5prime)) and canmerge_3prime)

    logging.debug("leaf exon=%s" % (str(leaf_exon)))
    logging.debug("match exon=%s" % (str(match_exon)))
    logging.debug("merge 5'=%s merge 3'=%s" % (merge_5prime, merge_3prime))

    if merge_5prime or merge_3prime:        
        merge_tuple = make_trimmed_merge_tuple(G, leaf_exon, leaf_score, 
                                               match_exon, match_score, 
                                               merge_5prime, merge_3prime)
        logging.debug("merge exon=%s" % (str(merge_tuple.node)))
        logging.debug("merge tuple=%s" % (str(merge_tuple)))
        
        if merge_tuple.node != leaf_exon:
            return merge_tuple, [leaf_exon] 
    return None, None


def merge_isoform_graph(G, mode="merge", overhang_threshold=15):
    # build interval tree
    exon_tree = IntervalTree()
    for node in G.nodes_iter():
        if node.is_exon:
            exon_tree.insert_interval(Interval(node.start, node.end, value=node))
    # get all unconnected nodes at start/end of graph
    start_nodes = [n for n,d in G.in_degree_iter()
                   if d == 0]
    end_nodes = [n for n,d in G.out_degree_iter()
                 if d == 0]
    leaf_exons = start_nodes + end_nodes

    merges_performed = 0
    merged_exon_data = []
    replaced_exons = set()
    # check each leaf exon for possible merges
    for leaf_exon in leaf_exons:
        # do a safety check to make sure algorithm is working
        if leaf_exon[3] == False:
            logging.critical("leaf is intron")
            sys.exit(0)        
        # intersect with all exons
        exon_hits = exon_tree.find(leaf_exon[0], leaf_exon[1])
        for exon_hit in exon_hits:
            match_exon = exon_hit.value
            # do not merge with self
            if leaf_exon == match_exon:
                continue
            # skip matches to incompatible strand
            if not cmp_strand(leaf_exon.strand, match_exon.strand):
                continue
            # perform specific types of merges according to mode
            merge_tuple = None
            if mode == "merge":
                merge_tuple, merged_exons = merge_intron_compatible_exons(G, leaf_exon, match_exon)
            elif mode == "trim":
                merge_tuple, merged_exons = merge_overhanging_exons(G, leaf_exon, match_exon, 
                                                                    overhang_threshold)
            if merge_tuple != None:
                #logging.debug("compare: %s %s" % (str(leaf_exon), str(match_exon)))
                #logging.debug("merge tuple: %s" % str(merge_tuple))
                #logging.debug("replaced exons: %s" % str(merged_exons))
                merged_exon_data.append(merge_tuple)
                merges_performed += 1
                replaced_exons.update(merged_exons)
    # remove all exons that were replaced by merges
    G.remove_nodes_from(replaced_exons)
    # reconstruct the graph by adding all the merged exons
    for merge_tuple in merged_exon_data:
        add_merged_exon(G, merge_tuple)
    logging.debug("merges performed=%d" % merges_performed)
    return merges_performed, G

def flatten_subgraph(G, nodes):
    cluster_tree = ClusterTree(0,1)
    strand = '.'
    id = 0
    for node in nodes:
        start, end, node_strand, is_exon = node
        if not is_exon:
            continue
        if strand == '.':
            strand = node_strand
        else:
            assert strand == node_strand
        cluster_tree.insert(node[0], node[1], id)
        id += 1
    flat_exons = []
    for region in cluster_tree.getregions():
        start, end, ids = region
        flat_exons.append((start,end))
    return strand, flat_exons

def flatten_locus(G):
    # break up the graph into connected components
    H = G.to_undirected()
    connected_nodes = nx.connected_components(H)
    for nodelist in connected_nodes:
        # one flattened transcript for each connected component
        yield flatten_subgraph(G, nodelist)

def find_all_paths(G, start_node, end_node, path=None):
    if path == None:
        path = []
    if start_node == end_node:
        yield path + [start_node]
    else:
        path.append(start_node)
        for next_node in G.successors(start_node):
            for newpath in find_all_paths(G, next_node, end_node, path):
                yield newpath
        path.pop()
    
def get_path_exons(path):
    strand = '.'
    exons = []
    for node in path:
        start, end, node_strand, is_exon = node
        if not is_exon:
            continue
        # figure out strand    
        if strand == '.':
            strand = node_strand
        else:
            assert strand == node_strand        
        exons.append((start, end))
    return strand, exons

def get_path_score(G, path):
    return sum([G.node[u][GRAPH_SCORE_KEY] for u in path])

def initialize_path_graph(G):
    # add dummy nodes to start and end
    dummystart = get_start_node()
    dummyend = get_end_node()
    # set edge weights
    for node in G.nodes():
        score = G.node[node][GRAPH_SCORE_KEY]
        # connect start nodes to dummy intron node
        if G.in_degree(node) == 0:
            G.add_edge(dummystart, node)
        if G.out_degree(node) == 0:
            G.add_edge(node, dummyend)
        # convert node weight to edge weights
        for pred in G.predecessors(node):
            G.edge[pred][node]['weight'] = score
        # outgoing edges go to introns, which have zero score
        for succ in G.successors(node):
            G.edge[node][succ]['weight'] = 0.0
    return dummystart, dummyend

def get_connected_components(G):
    H = G.to_undirected()
    connected_nodes = nx.connected_components(H)
    for nodelist in connected_nodes:
        yield nx.DiGraph(G.subgraph(nodelist))

def get_high_scoring_paths(G, start_node, end_node,
                           fraction_major_isoform,
                           max_isoforms,
                           max_paths,
                           clip_start=True,
                           clip_end=True):
    # get all paths in this subgraph
    paths = []
    for path in find_all_paths(G, start_node, end_node):
        # trim dummy start and end nodes
        if clip_start:
            path = path[1:]
        if clip_end:
            path = path[:-1]
        # get path score
        score = get_path_score(G, path)
        paths.append((score, path))
        if len(paths) > max_paths:
            logging.warning("Max paths >%d exceeded for this start/end node pair" % (max_paths))
            break
    # log
    logging.debug("%d paths" % len(paths))
    # rank paths
    paths = sorted(paths, key=operator.itemgetter(0), reverse=True)
    # choose the high scoring paths for each subgraph        
    major_isoform_score = paths[0][0]
    score_threshold = major_isoform_score * fraction_major_isoform
    # logging messages        
    logging.debug("major_isoform_score=%f score_threshold=%f" % 
                  (major_isoform_score, score_threshold))
    for i, path_tuple in enumerate(paths):
        logging.debug("rank=%d/%d score=%f threshold=%f major_isoform_score=%f" % 
                      (i, len(paths), path_tuple[0], score_threshold, major_isoform_score))
    paths = [path_tuple for path_tuple in paths
             if path_tuple[0] >= score_threshold]
    # log
    logging.debug("%d high scoring paths" % len(paths))
    # limit to max number of paths
    paths = paths[:max_isoforms]
    # yield paths
    for path_tuple in paths:
        score, path = path_tuple
        strand, exons = get_path_exons(path)
        yield strand, score, exons

    
def get_isoforms(G, fraction_major_isoform=0.50, 
                 max_isoforms=5,
                 max_paths=10000):  
    fraction_major_isoform = 0.99
    max_isoforms=1
    max_paths=10000
    # break up the graph into connected components
    for Gsub in get_connected_components(G):
        dummystart, dummyend = initialize_path_graph(Gsub)
        isoforms = list(get_high_scoring_paths(Gsub, dummystart, dummyend,
                                               fraction_major_isoform,
                                               max_isoforms,
                                               max_paths,
                                               clip_start=True,
                                               clip_end=True))
        yield isoforms

def get_tss_isoforms(G):
    fraction_major_isoform = 0.99
    max_isoforms=1
    max_paths=10000
    # break up the graph into connected components
    for Gsub in get_connected_components(G):
        start_nodes = [n for (n,d) in Gsub.in_degree_iter()
                       if (d == 0)]
        # convert node weights to edge weights
        # TODO: not using edge weights yet, but will in future
        dummystart, dummyend = initialize_path_graph(Gsub)
        tss_isoforms = []
        for start_node in start_nodes:
            tss_isoforms.extend(list(get_high_scoring_paths(Gsub, start_node, dummyend,
                                                            fraction_major_isoform,
                                                            max_isoforms,
                                                            max_paths,
                                                            clip_start=False,
                                                            clip_end=True)))
        yield tss_isoforms


def display_graph(G):
    # display graph
    pos = nx.spring_layout(G)
    plt.figure()
    nx.draw(G,pos)
    plt.show()

def merge_locus_transcripts(transcripts, mode="flatten"):
    '''
    modes:
    flatten = merge all isoforms together into a single transcript
    isoforms = output the top X isoforms
    tss = output one isoform (the highest ranked one) for each TSS
    '''
    # cluster transcripts to simplify merge
    transcripts = cluster_single_exon_transcripts(transcripts)
    # convert to graph form
    chrom = transcripts[0].chrom
    G = build_isoform_graph(transcripts)
    # iteratively merge the graph
    trim_rounds = 5
    merge_limit = 1e5
    for i in xrange(trim_rounds):
        merges_performed = 1
        while merges_performed > 0:
            merges_performed, G = merge_isoform_graph(G, mode="merge")
            if merges_performed > merge_limit:
                logging.debug("computational error could not merge transcripts")
                break
        trimmed, G = merge_isoform_graph(G, mode="trim")
        if not trimmed:
            break

    if mode == "flatten":
        merged_transcripts = []
        for strand, flat_exons in flatten_locus(G):
            strand_string = decode_strand(strand)
            merged_transcripts.append([Transcript(-1, chrom, strand_string, 0, flat_exons)])
    else:
        # make a transcript and use the path score
        merged_transcripts = []
        if mode == "isoforms":
            isoform_lists = list(get_isoforms(G))
        elif mode == "tss":
            isoform_lists = list(get_tss_isoforms(G))
        
        for gene_isoforms in isoform_lists:            
            isoform_transcripts = []      
            for strand, score, exons in gene_isoforms:
                strand_string = decode_strand(strand)
                t = Transcript(id=-1, 
                               chrom=chrom, 
                               strand=strand_string, 
                               score=score, 
                               exons=exons)
                isoform_transcripts.append(t)
            merged_transcripts.append(isoform_transcripts)        
    return merged_transcripts


def merge_transcripts(fileh, outfileh, mode):
    # counter used to assign unique ids to transcripts
    cur_tx_id = 0
    # parse by chromosome and locus
    for locus_features in parse_gtf_chrom_loci(fileh, cuffcompare.classify_gtf_attrs): 
        logging.debug("locus: %d transcripts" % (len(locus_features)))   
        # construct transcripts
        transcripts = [Transcript.from_gtf_features(features, score_attr="fpkm")
                       for features in locus_features]
        # logging (for debugging)
        for i, t in enumerate(transcripts):
            logging.debug("transcript %d: %s" % (i, t.to_bed12()))
        # merge transcripts in window together   
        merged_transcripts = merge_locus_transcripts(transcripts, mode)
        logging.debug("after merge: %d genes" % (len(merged_transcripts)))
        for gene_transcripts in merged_transcripts:
            for isoform_id, t in enumerate(gene_transcripts):
                t.id = 'TU_%07d_%d' % (cur_tx_id, isoform_id)
                print >>outfileh, t.to_bed12()
            cur_tx_id += 1

def get_modes():
    return ["flatten", "tss", "isoforms"]

if __name__ == '__main__':
    from optparse import OptionParser
    import sys
    logging.basicConfig(level=logging.DEBUG)
    parser = OptionParser()
    parser.add_option('--mode', dest="mode", default="flatten")
    options, args = parser.parse_args()    
    filename = args[0]
    merge_transcripts(open(filename), sys.stdout, options.mode)
