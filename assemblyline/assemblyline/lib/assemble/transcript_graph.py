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
from assemblyline.lib.base import GTFAttr

from base import NODE_SCORE, NODE_LENGTH, STRAND_SCORE, \
    TRANSCRIPT_IDS, SAMPLE_IDS, IS_REF
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

def get_transcript_node_map(G):
    t_node_map = collections.defaultdict(lambda: set())
    for n,d in G.nodes_iter(data=True):
        for t_id in d[TRANSCRIPT_IDS]:
            t_node_map[t_id].add(n)
    for t_id,nodes in t_node_map.iteritems():
        t_node_map[t_id] = sorted(nodes, key=operator.attrgetter('start'))
    return t_node_map

def add_node_undirected(G, n, t_id, sample_id, is_ref, strand, score):
    """
    add node to undirected graph
    
    each node in graph maintains attributes:
    'tids': set() of transcript id strings
    'sids': set() of sample id strings
    'ref': logical whether each strand is 'annotated'
    'strand_score': numpy array containing score on each strand
    'length': size of node in nucleotides
    """
    if n not in G: 
        attr_dict = {TRANSCRIPT_IDS: set(),
                     SAMPLE_IDS: set(),
                     IS_REF: np.zeros(3,bool),
                     NODE_LENGTH: (n.end - n.start),
                     STRAND_SCORE: np.zeros(3,float)} 
        G.add_node(n, attr_dict=attr_dict)
    nd = G.node[n]
    nd[TRANSCRIPT_IDS].add(t_id)
    if sample_id is not None:
        nd[SAMPLE_IDS].add(sample_id)
    nd[IS_REF][strand] |= is_ref
    nd[STRAND_SCORE][strand] += score

def create_undirected_transcript_graph(transcripts, gtf_sample_attr=None):
    '''
    add all transcripts to a single undirected graph
    '''
    # find the intron domains of the transcripts
    boundaries = find_exon_boundaries(transcripts)
    # initialize transcript graph as undirected at first
    G = nx.Graph()
    # add transcripts
    for t in transcripts:
        # get transcript attributes
        t_id = t.attrs[GTFAttr.TRANSCRIPT_ID]
        sample_id = t.attrs.get(gtf_sample_attr, None)
        is_ref = bool(int(t.attrs[GTFAttr.REF]))
        # split exons that cross boundaries and to get the
        # nodes in the transcript path
        nodes = []
        for exon in t.exons:
            for start,end in split_exon(exon, boundaries):
                nodes.append(Exon(start, end))
        # add nodes/edges to graph
        u = nodes[0]
        add_node_undirected(G, u, t_id, sample_id, is_ref, 
                            t.strand, t.score)
        for v in nodes[1:]:
            add_node_undirected(G, v, t_id, sample_id, is_ref, 
                                t.strand, t.score)
            G.add_edge(u, v)
            u = v
    return G

def redistribute_unstranded_transcripts(G, transcripts, transcript_node_map):
    """
    predicts strandedness of transcripts that lack strand information
    based on the distribution of fwd/rev strand transcripts across 
    overlapping nodes in the graph

    updates strand score attribute of transcript and strand score
    attribute of graph nodes
    """
    # iterate through transcripts and redistribute score
    num_resolved = 0
    unresolved = []
    node_score_delta_dict = collections.defaultdict(lambda: np.zeros(3,float))
    for t in transcripts:
        # ignore stranded transcripts
        if t.strand != NO_STRAND:
            continue
        nodes = transcript_node_map[t.attrs[GTFAttr.TRANSCRIPT_ID]]
        # sum the coverage score across the transcript
        score_arr = np.zeros(3,float)
        for n in nodes:
            score_arr += G.node[n][STRAND_SCORE]
        # calculate total mass across transcript
        total_strand_score = score_arr[POS_STRAND] + score_arr[NEG_STRAND]
        # if there is "stranded" mass on any of the nodes comprising
        # this transcript, then use the proportion of fwd/rev mass
        # to redistribute
        if total_strand_score > 0:
            # proportionally assign unstranded mass based on amount of
            # plus and minus strand mass
            pos_frac = score_arr[POS_STRAND] / float(total_strand_score)
            score_delta_arr = t.score * np.array((pos_frac, 1.0-pos_frac, -1.0))
            # save all score adjustments in a dictionary and wait to apply
            # until strand fractions are computed for all transcripts
            for n in nodes:
                node_score_delta_dict[n] += score_delta_arr
            # adjust strand score for transcript
            t.attrs[STRAND_SCORE] = t.score * np.array((pos_frac, 1.0-pos_frac, 0.0)) 
            num_resolved += 1
        else:
            unresolved.append(t)
    # subtract score from unstranded and add to stranded
    for n, score_delta_arr in node_score_delta_dict.iteritems():
        G.node[n][STRAND_SCORE] += score_delta_arr
    return num_resolved, unresolved

def redistribute_unstranded_node_clusters(G, transcripts, transcript_node_map):
    # find set of unresolved nodes
    unresolved_nodes = set()
    for t in transcripts:
        unresolved_nodes.update(transcript_node_map[t.attrs[GTFAttr.TRANSCRIPT_ID]])
    unresolved_nodes = sorted(unresolved_nodes, key=operator.attrgetter('start'))
    # cluster unresolved nodes
    cluster_tree = ClusterTree(0,1)
    for i,n in enumerate(unresolved_nodes):
        cluster_tree.insert(n.start, n.end, i)
    # try to assign stranded score to clusters
    node_score_delta_dict = collections.defaultdict(lambda: np.zeros(3,float))
    for start, end, indexes in cluster_tree.getregions():
        nodes = [unresolved_nodes[i] for i in indexes]
        # sum score across cluster
        score_arr = np.zeros(3,float)
        for n in nodes:
            score_arr += G.node[n][STRAND_SCORE]
        # calculate total mass across cluster
        total_strand_score = score_arr[POS_STRAND] + score_arr[NEG_STRAND]            
        if total_strand_score > 0:
            # proportionally assign unstranded mass based on amount of
            # plus and minus strand mass
            pos_frac = score_arr[POS_STRAND] / float(total_strand_score)
            score_delta_arr = t.score * np.array((pos_frac, 1.0-pos_frac, -1.0))
            # save all score adjustments in a dictionary and wait to apply
            # until strand fractions are computed for all transcripts
            for i,n in enumerate(nodes):
                node_score_delta_dict[n] += score_delta_arr
    del cluster_tree
    # subtract score from unstranded and add to stranded
    for n, score_delta_arr in node_score_delta_dict.iteritems():
        G.node[n][STRAND_SCORE] += score_delta_arr

def redistribute_score(G, transcripts):
    '''
    reallocate coverage mass of unstranded transcripts to the fwd/rev
    strand proportionately.
        
    after this method no nodes should contain both stranded and
    unstranded score
    '''
    # build a mapping from transcripts to graph nodes using the 
    # transcript id attributes of the nodes
    transcript_node_map = get_transcript_node_map(G)
    unresolved = []
    for t in transcripts:
        # set initial value in score vector
        strand_score = np.zeros(3,float)
        strand_score[t.strand] = t.score
        t.attrs[STRAND_SCORE] = strand_score
        if t.strand == NO_STRAND:
            unresolved.append(t)
    # try to reassign unstranded transcript score
    logging.debug("\t\tRedistributing score across unstranded transcripts")
    num_redist, unresolved = redistribute_unstranded_transcripts(G, unresolved, transcript_node_map)
    logging.debug("\t\tRescued %d unstranded transcripts (%d unresolved)" % (num_redist, len(unresolved)))
    if len(unresolved) > 0:
        # cluster remaining unstranded nodes and redistribute score
        # across the clusters
        logging.debug("\t\tRedistributing score across unstranded node clusters")
        redistribute_unstranded_node_clusters(G, unresolved, transcript_node_map)
        # try to reassign unstranded transcript score
        num_redist, unresolved = redistribute_unstranded_transcripts(G, unresolved, transcript_node_map)
        logging.debug("\t\tRescued another %d unstranded transcripts (%d unresolved)" % (num_redist, len(unresolved)))

def create_strand_transcript_maps(transcripts, gtf_sample_attr):
    """
    builds an undirected graph of all transcripts in order to
    resolve strandedness for unstranded transcripts and reallocate
    transcript score onto specific strands
    
    returns a list of dictionaries for each strand, dictionary
    keyed by transcript id
    """
    # build the initial transcript graph
    Gundir = create_undirected_transcript_graph(transcripts, gtf_sample_attr)
    # reallocate unstranded transcripts to fwd/rev strand according
    # fraction of fwd/rev score across transcript nodes
    redistribute_score(Gundir, transcripts)
    # now that score has been allocated, partition transcripts 
    # into fwd/rev/unknown strands
    transcript_maps = [{}, {}, {}]    
    for t in transcripts:
        for strand in xrange(0,3):
            if t.attrs[STRAND_SCORE][strand] > 1e-8:
                transcript_maps[strand][t.attrs[GTFAttr.TRANSCRIPT_ID]] = t
    # done with the undirected graph
    Gundir.clear()
    return transcript_maps

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

def add_transcript_directed(G, strand, boundaries, transcript):
    # get strand-specific score
    score = transcript.attrs[STRAND_SCORE][strand]
    # split exons that cross boundaries and to get the
    # nodes in the transcript path
    nodes = []
    for exon in transcript.exons:
        for start, end in split_exon(exon, boundaries):
            nodes.append(Exon(start,end))
    if strand == NEG_STRAND:
        nodes.reverse()
    # add nodes/edges to graph
    u = nodes[0]
    add_node_directed(G, u, transcript.attrs[GTFAttr.TRANSCRIPT_ID], score)
    for v in nodes[1:]:
        add_node_directed(G, v, transcript.attrs[GTFAttr.TRANSCRIPT_ID], score)
        G.add_edge(u,v)
        u = v
    
def create_strand_specific_graph(strand, transcripts):
    '''
    build strand-specific subgraphs for forward/reverse/unknown strands

    overhang_threshold: integer greater than zero specifying the 
    maximum exon overhang that can be trimmed to match an intron
    boundary.  exons that overhang more than this will be 
    considered independent transcript start sites or end sites
   
    nodes in the graph will have a 'score' attribute equal to the
    total strand-specific score for that node
    
    returns a 3-tuple containing DiGraph objects corresponding
    to the forward/reverse/unknown strand
    '''
    # find the intron domains of the transcripts
    boundaries = find_exon_boundaries(transcripts)
    # initialize transcript graph
    G = nx.DiGraph()
    # add transcripts
    for t in transcripts:
        add_transcript_directed(G, strand, boundaries, t)
    return G

def create_transcript_graphs(transcripts, gtf_sample_attr):
    '''
    nodes have the following attributes:
    chain: list of children nodes
    score: node weight
    length: total length of node

    generates (graph, strand, strand_transcript_map) tuples with transcript 
    graphs
    '''
    # redistribute transcript score by strand
    logging.debug("\tResolving unstranded transcripts")
    transcript_maps = create_strand_transcript_maps(transcripts, gtf_sample_attr)
    # create strand-specific graphs using redistributed score
    logging.debug("\tCreating transcript graphs")
    for strand, strand_transcript_map in enumerate(transcript_maps):
        # create strand specific transcript graph
        G = create_strand_specific_graph(strand, strand_transcript_map.values())
        yield G, strand, strand_transcript_map

def prune_transcript_graph(G, strand, strand_transcript_map,
                           min_trim_length=0, 
                           trim_utr_fraction=0.0,
                           trim_intron_fraction=0.0):
    '''
    overhang_threshold: integer greater than zero specifying the 
    maximum exon overhang that can be trimmed to match an intron
    boundary.  exons that overhang more than this will be 
    considered independent transcript start sites or end sites

    trim_utr_fraction: float specifying the fraction of the average UTR
    coverage below which the ends of the UTR will be trimmed

    trim_intron_fraction: float specifying the fraction of the average 
    intron coverage below which intronic nodes will be removed
    '''
    # trim utrs and intron retentions
    trim_graph(G, strand, min_trim_length, trim_utr_fraction, 
               trim_intron_fraction)
    # collapse consecutive nodes in graph
    H = collapse_strand_specific_graph(G)
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
            # adjust transcript score to account for changes in 
            # transcript length by renormalizing by length
            t = strand_transcript_map[t_id]
            score = t.attrs[STRAND_SCORE][strand]
            path_score_dict[tuple(nodes)] += score
            #path_length = sum(Gsub.node[n][NODE_LENGTH] for n in nodes)
            #new_score = score * min(1.0, float(t.length) / path_length)
            #path_score_dict[tuple(nodes)] += new_score
        yield Gsub, strand, path_score_dict.items()
