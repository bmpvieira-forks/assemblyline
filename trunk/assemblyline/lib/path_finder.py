'''
Created on Dec 4, 2010

@author: mkiyer
'''
import logging
import networkx as nx
import collections
import heapq

from bx.intervals.intersection import Interval, IntervalTree
from bx.intervals.cluster import ClusterTree

from cNode import cmp_strand, strand_str_to_int
from isoform_graph import POS_STRAND, NEG_STRAND, NO_STRAND, INTRON, EXON

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
    #logging.debug('CURRENT_PATH=%s' % str(path))
    if G.out_degree(start_node) == 0:
        yield path + [start_node]
        #logging.debug('END NODE REACHED=%s FINAL PATH=%s' % (str(start_node), str(path)))
    else:
        # search for partial paths starting from this node
        if start_node in partial_path_dict:
            for partial_path in partial_path_dict[start_node]:
                #logging.debug('Adding partial path=%s' % str(partial_path))                
                for newpath in find_paths_from_start(G, partial_path_dict, partial_path[-1], path + list(partial_path[:-1])):
                    yield newpath
        else:
            #logging.debug('ADDING NODE=%s' % str(start_node))
            # no partial path exists so use all possible paths
            for next_node in G.successors(start_node):
                #logging.debug('NEXT NODE=%s' % str(next_node))
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
    return 1.0e3 * total_score / total_length

def find_best_paths(G, partial_paths, id_score_map, num_paths=10, 
                    max_stored_paths=1000, 
                    max_iterations=1000000):
    # make a dictionary so it is possible to look up 
    # partial paths based on the start node
    partial_path_dict = collections.defaultdict(lambda: set())
    for p in partial_paths:
        for i in xrange(len(p)-1):
            partial_path_dict[p[i]].add(p[i:])
            #append(p[i:])
    partial_path_dict = dict(partial_path_dict)
    gene_id = 0
    # break up the graph into connected subgraphs
    for Gsub in get_connected_components(G):
        # find start nodes
        start_nodes = [n for (n,d) in Gsub.in_degree_iter()
                       if (d == 0)]
        tss_id = 0
        for start_node in start_nodes:
            paths = []
            iters = 0
            for path in find_paths_from_start(Gsub, partial_path_dict, start_node):
                score = score_path(G, id_score_map, path)
                if len(paths) < max_stored_paths:
                    # add path to graph
                    heapq.heappush(paths, (score, path))
                else:
                    # add path and discard the lowest scoring path
                    heapq.heappushpop(paths, (score, path))
                iters += 1
                if iters == max_iterations:
                    logging.warning("Path scoring reached max iterations before enumerating all paths")
                    break
            # discard low scoring paths from the heapq
            while len(paths) > num_paths:
                heapq.heappop(paths)           
            yield gene_id, tss_id, paths
            tss_id += 1
        gene_id += 1


def get_isoforms(G, transcripts, paths_per_tss=5):
    # build interval tree
    exon_tree = IntervalTree()
    for node in G.nodes_iter():
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
        for path in thread_transcript(G, exon_tree, strand, exons):
            # ignore paths of length 1 (useless for threading)
            if len(path) > 1:
                threaded_paths.add(tuple(path))
#                # TODO: check path for now
#                for i in xrange(len(path)-1):
#                    for j in xrange(i+1, len(path)):
#                        if nx.shortest_path(G, path[i], path[j]) == False:
#                            logging.error("PATH ERROR")
#                            assert False
#            for x in path:
#                assert x.node_type != INTRON
    logging.debug("threaded %s paths" % len(threaded_paths))
    for p in threaded_paths:
        logging.debug("threaded path=%s" % str([str(n) for n in p]))
        logging.debug("threaded node types=%s" % str([str(n.node_type) for n in p]))
        logging.debug("threaded exon sizes=%s" % str([str(n.end - n.start) for n in p]))
        logging.debug("threaded exon starts=%s" % str([str(n.start) for n in p]))
    logging.debug("/THREADING")
    # using transcript-threaded paths as a guide find full paths
    # through the graph
    logging.debug("ENUMERATING PATHS")
    for gene_id,tss_id,paths in find_best_paths(G, threaded_paths, id_score_map,
                                                paths_per_tss):
        yield gene_id, tss_id, paths
    logging.debug("/ENUMERATING PATHS")
    #for gene_paths in find_paths(self.G, threaded_paths):
    # find non-redundant set of paths
    # gene_paths = set(tuple(p) for p in gene_paths)
    # scored_gene_paths = [(p, score_path(self.G, id_score_map, p)) for p in gene_paths]
    # yield scored_gene_paths
