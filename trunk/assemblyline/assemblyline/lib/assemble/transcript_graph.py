'''
Created on Dec 17, 2011

@author: mkiyer
'''
import networkx as nx
import logging
import collections
import operator
import bisect
import numpy as np

from assemblyline.lib.bx.cluster import ClusterTree
from assemblyline.lib.transcript import Exon, POS_STRAND, NEG_STRAND, NO_STRAND
from base import NODE_DENSITY, NODE_LENGTH, STRAND_DENSITY, TRANSCRIPT_IDS, init_directed_node_attrs
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

def add_node_undirected(G, n, t_id, strand, density):
    """
    add node to undirected graph
    
    each node in graph maintains attributes:
    'ids': set() of transcript id strings
    'strand_density': numpy array containing density on each strand
    'length': size of node in nucleotides
    """
    if n not in G: 
        attr_dict = {TRANSCRIPT_IDS: set(),
                     NODE_LENGTH: (n.end - n.start),
                     STRAND_DENSITY: np.zeros(3,float)} 
        G.add_node(n, attr_dict=attr_dict)
    nd = G.node[n]
    nd[TRANSCRIPT_IDS].add(t_id)
    nd[STRAND_DENSITY][strand] += density

def create_undirected_transcript_graph(transcripts):
    '''
    add all transcripts to a single undirected graph
    '''
    # find the intron domains of the transcripts
    boundaries = find_exon_boundaries(transcripts)
    # initialize transcript graph as undirected at first
    G = nx.Graph()
    # add transcripts
    for t in transcripts:
        # split exons that cross boundaries and to get the
        # nodes in the transcript path
        nodes = []
        for exon in t.exons:
            for start,end in split_exon(exon, boundaries):
                nodes.append(Exon(start, end))
        # add nodes/edges to graph
        u = nodes[0]
        add_node_undirected(G, u, t.id, t.strand, t.fpkm)
        for v in nodes[1:]:
            add_node_undirected(G, v, t.id, t.strand, t.fpkm)
            G.add_edge(u, v)
            u = v
    return G

def redistribute_unstranded_transcripts(G, transcripts, transcript_node_map):
    """
    predicts strandedness of transcripts that lack strand information
    based on the distribution of fwd/rev strand transcripts across 
    overlapping nodes in the graph

    updates strand density attribute of transcript and strand density
    attribute of graph nodes
    """
    # iterate through transcripts and redistribute density
    num_resolved = 0
    unresolved = []
    node_density_delta_dict = collections.defaultdict(lambda: np.zeros(3,float))
    for t in transcripts:
        # ignore stranded transcripts
        if t.strand != NO_STRAND:
            continue
        nodes = transcript_node_map[t.id]
        # sum the coverage density across the transcript
        density_arr = np.zeros(3,float)
        for n in nodes:
            density_arr += G.node[n][STRAND_DENSITY]
        # calculate total mass across transcript
        total_strand_density = density_arr[POS_STRAND] + density_arr[NEG_STRAND]
        # if there is "stranded" mass on any of the nodes comprising
        # this transcript, then use the proportion of fwd/rev mass
        # to redistribute
        if total_strand_density > 0:
            # proportionally assign unstranded mass based on amount of
            # plus and minus strand mass
            pos_frac = density_arr[POS_STRAND] / float(total_strand_density)
            density_delta_arr = t.density * np.array((pos_frac, 1.0-pos_frac, -1.0))
            # save all density adjustments in a dictionary and wait to apply
            # until strand fractions are computed for all transcripts
            for n in nodes:
                node_density_delta_dict[n] += density_delta_arr
            # adjust strand density for transcript
            t.attrs[STRAND_DENSITY] = t.density * np.array((pos_frac, 1.0-pos_frac, 0.0)) 
            num_resolved += 1
        else:
            unresolved.append(t)
    # subtract density from unstranded and add to stranded
    for n, density_delta_arr in node_density_delta_dict.iteritems():
        G.node[n][STRAND_DENSITY] += density_delta_arr
    return num_resolved, unresolved

def redistribute_unstranded_node_clusters(G, transcripts, transcript_node_map):
    # find set of unresolved nodes
    unresolved_nodes = set()
    for t in transcripts:
        unresolved_nodes.update(transcript_node_map[t.id])
    unresolved_nodes = sorted(unresolved_nodes, key=operator.attrgetter('start'))
    # cluster unresolved nodes
    cluster_tree = ClusterTree(0,1)
    for i,n in enumerate(unresolved_nodes):
        cluster_tree.insert(n.start, n.end, i)
    # try to assign stranded density to clusters
    node_density_delta_dict = collections.defaultdict(lambda: np.zeros(3,float))
    for start, end, indexes in cluster_tree.getregions():
        nodes = [unresolved_nodes[i] for i in indexes]
        # sum density across cluster
        density_arr = np.zeros(3,float)
        for n in nodes:
            density_arr += G.node[n][STRAND_DENSITY]
        # calculate total mass across cluster
        total_strand_density = density_arr[POS_STRAND] + density_arr[NEG_STRAND]            
        if total_strand_density > 0:
            # proportionally assign unstranded mass based on amount of
            # plus and minus strand mass
            pos_frac = density_arr[POS_STRAND] / float(total_strand_density)
            density_delta_arr = t.density * np.array((pos_frac, 1.0-pos_frac, -1.0))
            # save all density adjustments in a dictionary and wait to apply
            # until strand fractions are computed for all transcripts
            for i,n in enumerate(nodes):
                node_density_delta_dict[n] += density_delta_arr
    del cluster_tree
    # subtract density from unstranded and add to stranded
    for n, density_delta_arr in node_density_delta_dict.iteritems():
        G.node[n][STRAND_DENSITY] += density_delta_arr

def redistribute_density(G, transcripts):
    '''
    reallocate coverage mass of unstranded transcripts to the fwd/rev
    strand proportionately.
        
    after this method no nodes should contain both stranded and
    unstranded density
    '''
    # build a mapping from transcripts to graph nodes using the 
    # transcript id attributes of the nodes
    transcript_node_map = get_transcript_node_map(G)
    unresolved = []
    for t in transcripts:
        # set initial value in density vector
        strand_density = np.zeros(3,float)
        strand_density[t.strand] = t.density
        t.attrs[STRAND_DENSITY] = strand_density
        if t.strand == NO_STRAND:
            unresolved.append(t)
    # try to reassign unstranded transcript density
    logging.debug("\t\tRedistributing density across unstranded transcripts")
    num_redist, unresolved = redistribute_unstranded_transcripts(G, unresolved, transcript_node_map)
    logging.debug("\t\tRescued %d unstranded transcripts (%d unresolved)" % (num_redist, len(unresolved)))
    if len(unresolved) > 0:
        # cluster remaining unstranded nodes and redistribute density
        # across the clusters
        logging.debug("\t\tRedistributing density across unstranded node clusters")
        redistribute_unstranded_node_clusters(G, unresolved, transcript_node_map)
        # try to reassign unstranded transcript density
        num_redist, unresolved = redistribute_unstranded_transcripts(G, unresolved, transcript_node_map)
        logging.debug("\t\tRescued another %d unstranded transcripts (%d unresolved)" % (num_redist, len(unresolved)))

def create_strand_transcript_maps(transcripts):
    """
    builds an undirected graph of all transcripts in order to
    resolve strandedness for unstranded transcripts and reallocate
    transcript density onto specific strands
    
    returns a list of dictionaries for each strand, dictionary
    keyed by transcript id
    """
    # build the initial transcript graph
    Gundir = create_undirected_transcript_graph(transcripts)
    # reallocate unstranded transcripts to fwd/rev strand according
    # fraction of fwd/rev density across transcript nodes
    redistribute_density(Gundir, transcripts)
    # now that density has been allocated, partition transcripts 
    # into fwd/rev/unknown strands
    transcript_maps = [{}, {}, {}]    
    for t in transcripts:
        for strand in xrange(0,3):
            if t.attrs[STRAND_DENSITY][strand] > 1e-8:
                transcript_maps[strand][t.id] = t
    # done with the undirected graph
    Gundir.clear()
    return transcript_maps

def add_node_directed(G, n, t_id, density):
    """
    add node to directed (strand-specific) graph
    """
    if n not in G: 
        G.add_node(n, attr_dict=init_directed_node_attrs(n))
    nd = G.node[n]
    nd[TRANSCRIPT_IDS].add(t_id)
    nd[NODE_DENSITY] += density

def add_transcript_directed(G, strand, boundaries, transcript):
    # get strand-specific density
    density = transcript.attrs[STRAND_DENSITY][strand]
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
    add_node_directed(G, u, transcript.id, density)
    for v in nodes[1:]:
        add_node_directed(G, v, transcript.id, density)
        G.add_edge(u,v)
        u = v
    
def create_strand_specific_graph(strand, transcripts):
    '''
    build strand-specific subgraphs for forward/reverse/unknown strands

    overhang_threshold: integer greater than zero specifying the 
    maximum exon overhang that can be trimmed to match an intron
    boundary.  exons that overhang more than this will be 
    considered independent transcript start sites or end sites
   
    nodes in the graph will have a 'density' attribute equal to the
    total strand-specific density for that node
    
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

def create_transcript_graphs(transcripts):
    '''
    min_length: minimum transcript length
    
    
    nodes have the following attributes:
    chain: list of children nodes
    density: coverage density at node 
    length: total length of node

    generates (graph, strand, strand_transcript_map) tuples with transcript 
    graphs
    '''
    # redistribute transcript density by strand
    logging.debug("\tResolving unstranded transcripts")
    transcript_maps = create_strand_transcript_maps(transcripts)
    # create strand-specific graphs using redistributed density
    logging.debug("\tCreating transcript graphs")
    for strand, strand_transcript_map in enumerate(transcript_maps):
        # create strand specific transcript graph
        G = create_strand_specific_graph(strand, strand_transcript_map.values())
        yield G, strand, strand_transcript_map

def prune_transcript_graph(G, strand, strand_transcript_map,
                           overhang_threshold=0, 
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
    trim_graph(G, strand, overhang_threshold, trim_utr_fraction, 
               trim_intron_fraction)
    # collapse consecutive nodes in graph
    #H = collapse_strand_specific_graph(G)
    H = G
    # get connected components of graph which represent independent genes
    # unconnected components are considered different genes
    Gsubs = nx.weakly_connected_component_subgraphs(H)
    for Gsub in Gsubs:
        # get partial path data supporting graph
        transcript_node_map = get_transcript_node_map(Gsub)
        path_density_dict = collections.defaultdict(lambda: 0)
        for t_id, nodes in transcript_node_map.iteritems():
            # reverse path for negative strand transcripts
            if strand == NEG_STRAND:
                nodes.reverse()
            # adjust transcript density to account for changes in 
            # transcript length by renormalizing by length
            t = strand_transcript_map[t_id]
            density = t.attrs[STRAND_DENSITY][strand]
            #path_density_dict[tuple(nodes)] += t.attrs[STRAND_DENSITY][strand]
            path_length = sum(Gsub.node[n][NODE_LENGTH] for n in nodes)
            new_density = density * max(1.0, float(t.length) / path_length)
            path_density_dict[tuple(nodes)] += new_density
        yield Gsub, strand, path_density_dict.items()
