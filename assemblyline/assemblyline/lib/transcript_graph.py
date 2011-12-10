'''
Created on Feb 13, 2011

@author: mkiyer
'''
import networkx as nx
import logging
import collections
import operator
import bisect
import numpy as np

from transcript import Exon, POS_STRAND, NEG_STRAND, NO_STRAND, cmp_strand
from assembler_base import NODE_DENSITY, NODE_LENGTH, STRAND_DENSITY, \
    TRANSCRIPT_IDS

def find_intron_starts_and_ends(transcripts):
    '''
    input: list of Transcript objects
    output: intron starts/ends by strand (list to index strand, list of starts)    
    
    list of lists of intron starts, list of lists of intron ends
    '''
    # TODO: combine 'find_exon_boundaries' with this
    intron_starts = [set(), set(), set()]
    intron_ends = [set(), set(), set()]
    for transcript in transcripts:
        exons = transcript.exons
        strand = transcript.strand
        if len(exons) == 1:
            continue
        if transcript.strand == NO_STRAND:
            logging.warning("Multi-exonic unstranded transcript detected!")
            continue
        intron_starts[strand].add(exons[0].end)
        for e in transcript.exons[1:-1]:
            intron_starts[strand].add(e.end)
            intron_ends[strand].add(e.start)
        intron_ends[strand].add(exons[-1].start)
    # combine pos/neg to get unstranded
    intron_starts[NO_STRAND] = intron_starts[POS_STRAND].union(intron_starts[NEG_STRAND])
    intron_ends[NO_STRAND] = intron_ends[POS_STRAND].union(intron_ends[NEG_STRAND])
    # sort lists
    for strand in (POS_STRAND, NEG_STRAND, NO_STRAND):
        intron_starts[strand] = sorted(intron_starts[strand])
        intron_ends[strand] = sorted(intron_ends[strand])
    return intron_starts, intron_ends

def trim_left(start, end, intron_starts, intron_ends, overhang_threshold):
    trim_start = start
    # search for the nearest intron end greater than
    # or equal to the current exon start
    i = bisect.bisect_left(intron_ends, start)
    if i != len(intron_ends):
        nearest_intron_end = intron_ends[i]
        # cannot trim past end of exon and cannot
        # trim more than the overhang threshold
        if ((nearest_intron_end < end) and
            (nearest_intron_end <= start + overhang_threshold)):
            trim_start = nearest_intron_end
    return trim_start

def trim_right(start, end, intron_starts, intron_ends, overhang_threshold):
    trim_end = end
    if trim_right:
        # search for nearest intron start less than
        # or equal to the current exon end
        i = bisect.bisect_right(intron_starts, end)
        if i > 0:
            nearest_intron_start = intron_starts[i-1]
            # cannot trim past start of exon and cannot
            # trim more than the overhang threshold
            if ((nearest_intron_start > start) and
                (nearest_intron_start >= end - overhang_threshold)):
                trim_end = nearest_intron_start    
    return trim_end

def trim_transcript(exons, intron_starts, intron_ends, 
                    overhang_threshold):
    # only the very first and last exon of each transcript can be 
    # trimmed
    trim_start = trim_left(exons[0].start, exons[0].end, 
                           intron_starts, intron_ends, 
                           overhang_threshold)
    trim_end = trim_right(exons[-1].start, exons[-1].end, 
                          intron_starts, intron_ends, 
                          overhang_threshold)
    # if both start and end were trimmed it is
    # possible that they could be trimmed to 
    # match different introns and generate a 
    # weird and useless trimmed exon that is not
    # compatible with either the left or right
    # introns. resolve these situations by choosing
    # the smaller of the two trimmed distances as
    # the more likely
    if trim_end <= trim_start:
        left_trim_dist = trim_start - exons[0].start
        right_trim_dist = exons[-1].end - trim_end
        assert left_trim_dist >= 0
        assert right_trim_dist >= 0
        if left_trim_dist <= right_trim_dist:
            trim_start, trim_end = trim_start, exons[-1].end
        else:
            trim_start, trim_end = exons[0].start, trim_end
    assert trim_start < trim_end
    return trim_start, trim_end

def trim_transcripts(transcripts, overhang_threshold):
    if overhang_threshold <= 0:
        return
    # find the set of intron boundaries that govern 
    # trimming
    intron_starts, intron_ends = find_intron_starts_and_ends(transcripts)    
    for t in transcripts:
        trim_start, trim_end = trim_transcript(t.exons, 
                                               intron_starts[t.strand],
                                               intron_ends[t.strand],
                                               overhang_threshold)
        assert trim_start < trim_end
        # finally, modify the transcript
        t.exons[0].start = trim_start
        t.exons[-1].end = trim_end 

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

def add_node_to_graph(G, n, id, strand, density):
    """
    add node to graph
    
    each node in graph maintains attributes:
    'ids': set() of transcript id strings
    'strand_density': numpy array containing density on each strand
    'length': size of node in nucleotides
    """
    if n not in G: 
        attr_dict={TRANSCRIPT_IDS: set(),
                   NODE_LENGTH: (n.end - n.start),
                   STRAND_DENSITY: np.zeros(3,float)} 
        G.add_node(n, attr_dict=attr_dict)
    nd = G.node[n]
    nd[TRANSCRIPT_IDS].add(id)
    nd[STRAND_DENSITY][strand] += density

def add_edge_to_graph(G, u, v):
    """
    add edge to graph
    """
    if not G.has_edge(u, v):
        G.add_edge(u, v, attr_dict={})

def add_transcripts_to_graph(transcripts, overhang_threshold=0):
    '''
    overhang_threshold: integer greater than zero specifying the 
    maximum exon overhang that can be trimmed to match an intron
    boundary.  exons that overhang more than this will be 
    considered independent transcript start sites or end sites
    '''
    # trim the transcripts (modifies transcripts in place)
    trim_transcripts(transcripts, overhang_threshold)
    # find the intron domains of the transcripts
    boundaries = find_exon_boundaries(transcripts)
    # initialize transcript graph
    G = nx.DiGraph()
    # add transcripts
    for t in transcripts:
        if t.strand == NO_STRAND:
            # ensure unstranded transcripts have only one exon
            assert len(t.exons) == 1
        # split exons that cross boundaries and to get the
        # nodes in the transcript path
        nodes = []
        for exon in t.exons:
            nodes.extend([Exon(start, end) for start, end in 
                          split_exon(exon, boundaries)])
        # add nodes/edges to graph
        u = nodes[0]
        add_node_to_graph(G, u, t.id, t.strand, t.fpkm)
        for v in nodes[1:]:
            add_node_to_graph(G, v, t.id, t.strand, t.fpkm)
            # add edges if transcript is compatible
            if cmp_strand(t.strand, POS_STRAND):
                add_edge_to_graph(G, u, v)
            if cmp_strand(t.strand, NEG_STRAND):
                add_edge_to_graph(G, v, u)
            u = v
    return G

def get_transcript_node_map(G):
    t_node_map = collections.defaultdict(lambda: set())
    for n,d in G.nodes_iter(data=True):
        for id in d[TRANSCRIPT_IDS]:
            t_node_map[id].add(n)
    for id,nodes in t_node_map.iteritems():
        nodes = sorted(nodes, key=operator.attrgetter('start'))
        t_node_map[id] = nodes
    return t_node_map

def redistribute_iteration(G, transcripts):
    # build a mapping from transcripts to graph nodes using the 
    # transcript id attributes of the nodes
    transcript_node_map = get_transcript_node_map(G)
    # iterate through transcripts and redistribute density
    node_density_delta_dict = collections.defaultdict(lambda: np.zeros(3,float))
    unresolved = []
    num_resolved = 0
    for t in transcripts:
        # only process unstranded transcripts
        if t.strand != NO_STRAND:
            continue
        nodes = transcript_node_map[t.id]
        # sum the coverage mass across the transcript
        mass_arr = np.zeros(3,float)
        total_node_length = 0        
        for n in nodes:
            d = G.node[n]
            mass_arr += d[NODE_LENGTH] * d[STRAND_DENSITY]
            total_node_length += d[NODE_LENGTH]
        # calculate effective transcript density across new node length
        # might not match transcript length due to chaining of nodes
        assert t.length == total_node_length        
        #t_adj_density = t.density * t.length / float(total_node_length)
        # calculate total mass across transcript
        total_strand_mass = mass_arr[POS_STRAND] + mass_arr[NEG_STRAND]
        # if there is "stranded" mass on any of the nodes comprising
        # this transcript, then use the proportion of fwd/rev mass
        # to redistribute the unstranded mass 
        if total_strand_mass > 0:
            # proportionally assign unstranded mass based on amount of
            # plus and minus strand mass
            pos_frac = mass_arr[POS_STRAND] / float(total_strand_mass)
            density_delta_arr = t.density * np.array((pos_frac, 1.0-pos_frac, -1.0))
            for n in nodes:
                # save all density adjustments in a dictionary and wait to apply
                # until strand fractions are computed for all transcripts
                node_density_delta_dict[n] += density_delta_arr
            num_resolved += 1
        else:
            unresolved.append(t)
    # enact all the adjustments at once
    for n, density_delta_arr in node_density_delta_dict.iteritems():
        # subtract density from unstranded and add to stranded
        G.node[n][STRAND_DENSITY] += density_delta_arr
    return num_resolved, unresolved

def redistribute_unstranded_density(G, transcripts):
    '''
    reallocate coverage mass of unstranded transcripts to the fwd/rev
    strand proportionately.
        
    after this method no nodes should contain both stranded and
    unstranded density
    '''
    num_redist = 1
    remaining = transcripts
    while (num_redist > 0):
        num_redist, remaining = redistribute_iteration(G, remaining)

def create_strand_specific_graphs(G):
    '''
    build separate subgraphs of G for forward/reverse/unknown strands
   
    nodes in the graph will have a 'density' attribute equal to the
    total strand-specific density for that node

    edges (u,v) in the graph will have a 'density' attribute corresponding
    to the fraction of total coverage flowing out of node 'u' that goes 
    to 'v'
    
    NOTE: the unknown strand graph will not have any edges, just nodes
    
    returns a 3-tuple containing DiGraph objects corresponding
    to the forward/reverse/unknown strand
    '''
    # make graphs
    GG = (nx.DiGraph(), nx.DiGraph(), nx.DiGraph())
    # add nodes
    for n,d in G.nodes_iter(data=True):
        strand_density = d[STRAND_DENSITY]
        transcript_ids = d[TRANSCRIPT_IDS]
        for strand in (POS_STRAND, NEG_STRAND, NO_STRAND):
            density = strand_density[strand]
            # use a slightly non-zero value to overcome any numerical precision issues
            if density > 1e-8:
                GG[strand].add_node(n, attr_dict={NODE_DENSITY: density,
                                                  TRANSCRIPT_IDS: transcript_ids})
    # add edges
    for u,v,d in G.edges_iter(data=True):            
        # since the edge is directional from u -> v, we can 
        # now determine strand based on coordinates of nodes
        strand = int(u.start >= v.end)
        # get density of both nodes
        u_strand_density = G.node[u][STRAND_DENSITY]
        v_strand_density = G.node[v][STRAND_DENSITY]
        # check stranded density
        if ((u_strand_density[strand] > 1e-8) and 
            (v_strand_density[strand] > 1e-8)):
            GG[strand].add_edge(u, v, attr_dict=d)
        # ignore neg strand edges for graph on unknown strand
        # by convention (prevents bidirectional edges)
        if ((strand == POS_STRAND) and
            (u_strand_density[NO_STRAND] > 1e-8) and 
            (v_strand_density[NO_STRAND] > 1e-8)):
            GG[NO_STRAND].add_edge(u, v, attr_dict=d)
    return GG
    
def create_transcript_graph(transcripts, overhang_threshold=0):
    '''
    overhang_threshold: integer greater than zero specifying the 
    maximum exon overhang that can be trimmed to match an intron
    boundary.  exons that overhang more than this will be 
    considered independent transcript start sites or end sites
    
    returns a 3-tuple containing fwd/rev/unknown transcript graphs
    
    nodes have the following attributes:
    chain: list of children nodes
    density: coverage density at node 
    length: total length of node
    
    edges have the following attributes:
    density: coverage density flowing through edge
    '''
    # build the initial transcript graph
    G = add_transcripts_to_graph(transcripts, overhang_threshold)
    # reallocate unstranded transcripts to fwd/rev strand according
    # fraction of fwd/rev density across transcript nodes
    redistribute_unstranded_density(G, transcripts)
    # separate graph into strand-specific subgraphs
    HH = create_strand_specific_graphs(G)
    return HH