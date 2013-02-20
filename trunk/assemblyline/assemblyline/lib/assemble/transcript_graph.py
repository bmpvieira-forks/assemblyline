'''
Created on Dec 17, 2011

@author: mkiyer

AssemblyLine: transcriptome meta-assembly from RNA-Seq

Copyright (C) 2012 Matthew Iyer

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
'''
import networkx as nx
import logging
import collections
import operator
import bisect
import numpy as np

from assemblyline.lib.bx.cluster import ClusterTree
from assemblyline.lib.transcript import Exon, POS_STRAND, NEG_STRAND, NO_STRAND
from assemblyline.lib.base import GTFAttr, FLOAT_PRECISION
from base import NODE_SCORE, NODE_LENGTH, STRAND_SCORE, TRANSCRIPT_IDS 
from trim import trim_graph
from collapse import collapse_strand_specific_graph

def find_exon_boundaries(transcripts):
    '''
    input: a list of transcripts (not Node objects, these are transcripts)
    parsed directly from an input file and not added to an isoform graph
    yet. 
    
    output: sorted list of exon boundaries    
    '''
    exon_boundaries = set()
    # first add introns to the graph and keep track of
    # all intron boundaries
    for transcript in transcripts:
        # add transcript exon boundaries
        for exon in transcript.exons:
            # keep track of positions where introns can be joined to exons
            exon_boundaries.add(exon.start)
            exon_boundaries.add(exon.end)
    # sort the intron boundary positions and add them to interval trees
    return sorted(exon_boundaries)

def split_exon(exon, boundaries):
    """
    partition the exon given list of node boundaries
    
    generator yields (start,end) intervals for exon
    """
    # find the indexes into the intron boundaries list that
    # border the exon.  all the indexes in between these two
    # are overlapping the exon and we must use them to break
    # the exon into pieces 
    start_ind = bisect.bisect_right(boundaries, exon.start)
    end_ind = bisect.bisect_left(boundaries, exon.end)
    exon_splits = [exon.start] + boundaries[start_ind:end_ind] + [exon.end]
    for j in xrange(1, len(exon_splits)):
        start, end = exon_splits[j-1], exon_splits[j]
        yield start, end

def split_exons(t, boundaries):
    # split exons that cross boundaries and to get the
    # nodes in the transcript path
    nodes = []
    for exon in t.exons:
        for start,end in split_exon(exon, boundaries):
            nodes.append(Exon(start, end))
    return nodes

def resolve_strand(nodes, node_data):
    # find strand with highest score or strand
    # best supported by reference transcripts
    total_scores = [0.0, 0.0]
    ref_bp = [0, 0]
    for n in nodes:
        length = (n.end - n.start)
        nd = node_data[n]
        scores = nd['scores']
        total_scores[POS_STRAND] += (scores[POS_STRAND] * length)
        total_scores[NEG_STRAND] += (scores[NEG_STRAND] * length)
        ref_strands = nd['ref_strands']
        if ref_strands[POS_STRAND]:
            ref_bp[POS_STRAND] += length
        if ref_strands[NEG_STRAND]:
            ref_bp[NEG_STRAND] += length
    if sum(total_scores) > FLOAT_PRECISION:
        if total_scores[POS_STRAND] >= total_scores[NEG_STRAND]:
            return POS_STRAND
        else:
            return NEG_STRAND
    if sum(ref_bp) > 0:
        if ref_bp[POS_STRAND] >= ref_bp[NEG_STRAND]:
            return POS_STRAND
        else:
            return NEG_STRAND
    return NO_STRAND

def add_transcript(t, nodes, strand_transcript_maps, node_data):
    for n in nodes:
        nd = node_data[n]
        nd['scores'][t.strand] += t.score
    t_id = t.attrs[GTFAttr.TRANSCRIPT_ID]
    strand_transcript_maps[t.strand][t_id] = t

def partition_transcripts_by_strand(transcripts):
    # divide transcripts into independent regions of
    # transcription with a single entry and exit point
    boundaries = find_exon_boundaries(transcripts)
    node_data_func = lambda: {'ref_strands': [False, False],
                              'scores': [0.0, 0.0, 0.0]}
    node_data = collections.defaultdict(node_data_func)
    stranded_transcript_maps = [{}, {}, {}]
    stranded_ref_transcripts = [[], []]
    unresolved_transcripts = []
    for t in transcripts:
        # split exons that cross boundaries and get the
        # nodes that made up the transcript
        nodes = split_exons(t, boundaries)        
        is_ref = bool(int(t.attrs.get(GTFAttr.REF, "0")))
        if is_ref:
            # label nodes by ref strand
            for n in nodes:
                nd = node_data[n]
                nd['ref_strands'][t.strand] = True
            stranded_ref_transcripts[t.strand].append(t)
        elif t.strand != NO_STRAND:
            add_transcript(t, nodes, stranded_transcript_maps, node_data)
        else:
            unresolved_transcripts.append((t,nodes))
    # resolve unstranded transcripts
    logging.debug("\t\t%d unstranded transcripts" % 
                  (len(unresolved_transcripts)))
    if len(unresolved_transcripts) > 0:
        resolved = []
        still_unresolved_transcripts = []
        for t,nodes in unresolved_transcripts:
            t.strand = resolve_strand(nodes, node_data)
            if t.strand != NO_STRAND:
                resolved.append((t,nodes))
            else:
                still_unresolved_transcripts.append((t,nodes))
        for t,nodes in resolved:
            add_transcript(t, nodes, stranded_transcript_maps, node_data)
        unresolved_transcripts = still_unresolved_transcripts
    if len(unresolved_transcripts) > 0:
        logging.debug("\t\t%d unresolved transcripts" % 
                      (len(unresolved_transcripts)))
        # if there are still unresolved transcripts then we can try to
        # extrapolate and assign strand to clusters of nodes at once, as
        # long as some of the nodes have a strand assigned
        # cluster unresolved nodes
        unresolved_nodes = set()
        for t,nodes in unresolved_transcripts:
            unresolved_nodes.update(nodes)
        unresolved_nodes = sorted(unresolved_nodes, key=operator.attrgetter('start'))
        cluster_tree = ClusterTree(0,1)
        for i,n in enumerate(unresolved_nodes):
            cluster_tree.insert(n.start, n.end, i)
        # try to assign strand to clusters of nodes
        node_strand_map = {}
        for start, end, indexes in cluster_tree.getregions():
            nodes = [unresolved_nodes[i] for i in indexes]
            strand = resolve_strand(nodes, node_data)
            for n in nodes:
                node_strand_map[n] = strand
        # for each transcript assign strand to the cluster with
        # the best overlap
        unresolved_count = 0
        for t,nodes in unresolved_transcripts:
            strand_bp = [0,0]
            for n in nodes:
                strand = node_strand_map[n]
                if strand != NO_STRAND:
                    strand_bp[strand] += (n.end - n.start)
            total_strand_bp = sum(strand_bp)
            if total_strand_bp > 0:
                if strand_bp[POS_STRAND] >= strand_bp[NEG_STRAND]:
                    t.strand = POS_STRAND
                else:
                    t.strand = NEG_STRAND
            else:
                unresolved_count += 1
            add_transcript(t, nodes, stranded_transcript_maps, node_data)
        logging.debug("\t\tCould not resolve %d transcripts" % 
                      (unresolved_count))
        del cluster_tree    
    return stranded_transcript_maps, stranded_ref_transcripts

def add_node_directed(G, n, t_id, score):
    """
    add node to directed (strand-specific) graph
    """
    if n not in G: 
        G.add_node(n, attr_dict={TRANSCRIPT_IDS: set(),
                                 NODE_LENGTH: (n.end - n.start), 
                                 NODE_SCORE: 0.0})
    nd = G.node[n]
    nd[TRANSCRIPT_IDS].add(t_id)
    nd[NODE_SCORE] += score

def create_directed_graph(strand, transcripts):
    '''
    build strand-specific graph
    '''
    # find the intron domains of the transcripts
    boundaries = find_exon_boundaries(transcripts)
    # initialize transcript graph
    G = nx.DiGraph()
    # add transcripts
    for t in transcripts:
        t_id = t.attrs[GTFAttr.TRANSCRIPT_ID]
        # split exons that cross boundaries and get the
        # nodes that made up the transcript
        nodes = split_exons(t, boundaries)
        if strand == NEG_STRAND:
            nodes.reverse()
        # add nodes/edges to graph
        u = nodes[0]
        add_node_directed(G, u, t_id, t.score)
        for v in nodes[1:]:
            add_node_directed(G, v, t_id, t.score)
            G.add_edge(u,v)
            u = v
    return G

def create_transcript_graphs(transcripts):
    '''
    generates (graph, strand, transcript_map) tuples with transcript 
    graphs
    '''
    # partition transcripts by strand and resolve unstranded transcripts
    logging.debug("\tResolving unstranded transcripts")
    strand_transcript_maps, strand_ref_transcripts = \
        partition_transcripts_by_strand(transcripts)
    # create strand-specific graphs using redistributed score
    logging.debug("\tCreating transcript graphs")
    for strand, transcript_map in enumerate(strand_transcript_maps):
        # create strand specific transcript graph
        G = create_directed_graph(strand, transcript_map.values())
        yield G, strand, transcript_map

def get_transcript_node_map(G, node_attr=TRANSCRIPT_IDS):
    t_node_map = collections.defaultdict(lambda: set())
    for n,nd in G.nodes_iter(data=True):
        for t_id in nd[node_attr]:
            t_node_map[t_id].add(n)
    for t_id,nodes in t_node_map.iteritems():
        t_node_map[t_id] = sorted(nodes, key=operator.attrgetter('start'))
    return t_node_map

def prune_transcript_graph(G, strand, transcript_map,
                           min_trim_length=0, 
                           trim_utr_fraction=0.0,
                           trim_intron_fraction=0.0):
    '''
    trim_utr_fraction: float specifying the fraction of the average UTR
    coverage below which the ends of the UTR will be trimmed

    trim_intron_fraction: float specifying the fraction of the average 
    intron coverage below which intronic nodes will be removed
    '''
    # trim utrs and intron retentions
    trim_nodes = trim_graph(G, strand, min_trim_length, 
                            trim_utr_fraction, 
                            trim_intron_fraction)
    G.remove_nodes_from(trim_nodes)
    # collapse consecutive nodes in graph
    H = collapse_strand_specific_graph(G, transcript_map)
    # get connected components of graph which represent independent genes
    # unconnected components are considered different genes
    Gsubs = nx.weakly_connected_component_subgraphs(H)
    for Gsub in Gsubs:
        # get partial path data supporting graph
        transcript_node_map = get_transcript_node_map(Gsub)
        path_score_dict = collections.defaultdict(lambda: 0)
        for t_id, nodes in transcript_node_map.iteritems():
            # reverse path for negative strand transcripts
            if strand == NEG_STRAND:
                nodes.reverse()
            # get transcript scores
            t = transcript_map[t_id]
            path_score_dict[tuple(nodes)] += t.score
        yield Gsub, strand, path_score_dict.items()
