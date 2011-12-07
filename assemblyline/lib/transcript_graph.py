'''
Created on Feb 13, 2011

@author: mkiyer
'''
import networkx as nx
import logging
import collections
import itertools
import operator
import bisect
import numpy as np

from bx.intersection import Interval, IntervalTree
from bx.cluster import ClusterTree
from transcript import Exon, POS_STRAND, NEG_STRAND, NO_STRAND, \
    cmp_strand, strand_int_to_str

class TranscriptData(object):
    __slots__ = ('id', 'strand', 'density')
    def __init__(self, id, strand, density):
        self.id = id
        self.strand = strand
        self.density = float(density)
    def __repr__(self):
        return ("<%s(id='%s',strand='%s',density='%s')>" % 
                (self.__class__.__name__, self.id, self.strand, 
                 str(self.density))) 
    def __str__(self):
        return ("[id=%s,strand=%s,density=%.6f]" % 
                (self.id, strand_int_to_str(self.strand), self.density))

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
    #logging.debug('strand=%d txstart=%d txend=%d trimstart=%d trimend=%d' % 
    #              (t.strand, t.tx_start, t.tx_end, trim_start, trim_end))
    # if both start and end were trimmed it is
    # possible that they could be trimmed to 
    # match different introns and generate a 
    # weird and useless trimmed exon that is not
    # compatible with either the left or right
    # introns. resolve these situations by choosing
    # the smaller of the two trimmed distances as
    # the more likely
    if trim_end <= trim_start:
        #logging.warning('Trimming produced txstart=%d txend=%d trimstart=%d trimend=%d' % 
        #                (exons[0].start, exons[-1].end, trim_start, trim_end))
        left_trim_dist = trim_start - exons[0].start
        right_trim_dist = exons[-1].end - trim_end
        assert left_trim_dist >= 0
        assert right_trim_dist >= 0
        if left_trim_dist <= right_trim_dist:
            trim_start, trim_end = trim_start, exons[-1].end
        else:
            trim_start, trim_end = exons[0].start, trim_end
        #logging.warning('\tresolved to trimstart=%d trimend=%d' % (trim_start, trim_end))
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

def get_exon_node_edges(exon, tdata, boundaries):
    """
    exon: Exon object to split
    tdata: TranscriptData object
    
    returns lists of (node,data) and (edge,data) tuples
    """
    nodes = []
    edges = []
    for start, end in split_exon(exon, boundaries):
        # make node
        n = Exon(start, end)
        exon_tdata = TranscriptData(id=tdata.id, 
                                    strand=tdata.strand, 
                                    density=tdata.density)
        # make edges between split exon according to 
        # strand being assembled.
        if len(nodes) > 0:
            if cmp_strand(tdata.strand, NEG_STRAND):
                edges.append((n, nodes[-1][0], tdata))
            if cmp_strand(tdata.strand, POS_STRAND):
                edges.append((nodes[-1][0], n, tdata))
        nodes.append((n, exon_tdata))
    assert nodes[-1][0].end == exon.end
    return nodes, edges

def get_transcript_nodes_edges(transcript, boundaries):
    """
    transcript: Transcript object
    boundaries: list of exon/intron boundaries
    
    returns list of (node,data) tuples and list of (u,v,data) edge tuples
    """
    exons = transcript.exons
    strand = transcript.strand
    tdata = TranscriptData(id=transcript.id, strand=strand, 
                           density=transcript.fpkm)
    # split exons that cross boundaries and get the
    # new exons of the path
    tnodes = []
    tedges = []
    # add the first exon to initialize the loop
    # (all transcripts must have at least one exon)
    enodes, eedges = get_exon_node_edges(exons[0], tdata, boundaries)
    tnodes.extend(enodes)
    tedges.extend(eedges)
    for e2 in exons[1:]:
        # add exon
        enodes, eedges = get_exon_node_edges(e2, tdata, boundaries)
        # add intron -> exon edges
        if strand != NO_STRAND:
            if strand == NEG_STRAND:
                tedges.append((enodes[0][0], tnodes[-1][0], tdata))
            else:
                tedges.append((tnodes[-1][0], enodes[0][0], tdata))
        # update nodes/edges
        tnodes.extend(enodes)
        tedges.extend(eedges)
    return tnodes, tedges

def add_node_to_graph(G, n, tdata):
    if n not in G:  
        G.add_node(n, data={})        
    nd = G.node[n]
    nd['data'][tdata.id] = tdata

def add_edge_to_graph(G, u, v, tdata):
    if not G.has_edge(u, v):
        G.add_edge(u, v, data={})
    ed = G.edge[u][v]
    ed['data'][tdata.id] = tdata

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
    # keep mapping of transcripts to graph nodes
    transcript_graph_map = {}
    unstranded_ids = []
    # initialize transcript graph
    G = nx.DiGraph()
    # add transcripts
    for t in transcripts:
        if t.strand == NO_STRAND:
            # ensure unstranded transcripts have only one exon
            assert len(t.exons) == 1
            # save NO_STRAND transcripts for reprocessing later
            unstranded_ids.append(t.id)
        # add transcripts to graph
        node_data_list, edge_data_list = get_transcript_nodes_edges(t, boundaries)
        transcript_graph_map[t.id] = (node_data_list, edge_data_list)
        for n,nd in node_data_list:
            add_node_to_graph(G, n, nd)        
        for u,v,ed in edge_data_list:
            add_edge_to_graph(G, u, v, ed)
    return G, transcript_graph_map, unstranded_ids

def sum_strand_densities(G):
    '''
    compute the total density per strand at each node/edge

    a 'strand_density' numpy array will exist at each node/edge in the 
    graph with this information    
    '''
    for n,d in G.nodes_iter(data=True):
        # access list of tdata objects
        tdata_list = d['data'].values()
        # sum density at each node
        strand_density = np.zeros(3, float)
        for tdata in tdata_list:
            strand_density[tdata.strand] += tdata.density
        d['strand_density'] = strand_density
    for u,v,d in G.edges_iter(data=True):            
        # access list of tdata objects
        tdata_list = d['data'].values()
        strand_density = np.zeros(3, float)
        for tdata in tdata_list:
            strand_density[tdata.strand] += tdata.density
        d['strand_density'] = strand_density

def redistribute_unstranded_density(G, transcript_ids, transcript_graph_map):
    '''
    reallocate coverage mass of unstranded transcripts to the fwd/rev
    strand proportionately.  must call sum_strand_densities(G) to set
    the 'strand_density' attribute at each node/edge before calling this
    method
        
    after this method no nodes should contain both stranded and
    unstranded density
    '''
    # compute the fraction of fwd/rev coverage for each transcript
    for id in transcript_ids:
        tnodes,tedges = transcript_graph_map[id]        
        # sum the coverage mass across the transcript
        strand_mass = np.zeros(3,float)
        for n,d in tnodes:
            # TODO: can remove this, but check strand and density            
            assert d.strand == NO_STRAND
            assert d.density == G.node[n]['data'][id].density
            strand_mass += (n.end - n.start) * G.node[n]['strand_density']            
        # calculate total mass across transcript
        total_strand_mass = strand_mass[POS_STRAND] + strand_mass[NEG_STRAND]        
        # if there is "stranded" mass on any of the nodes comprising
        # this transcript, then use the proportion of fwd/rev mass
        # to re-apportion the unstranded mass 
        if total_strand_mass > 0:
            # proportionally assign unstranded mass based on amount of
            # plus and minus strand mass
            pos_frac = strand_mass[POS_STRAND] / float(total_strand_mass)
            redistribute_arr = np.array((pos_frac, 1.0 - pos_frac, -1.0))
            #density_adjustment = t.density * np.array((pos_frac, 1.0 - pos_frac, -1.0), dtype=float)           
            for n,d in tnodes:
                # subtract density from unstranded and add to stranded
                G.node[n]['strand_density'] += d.density * redistribute_arr
            for u,v,d in tedges:
                # TODO: can remove this, but check strand and density
                assert d.strand == NO_STRAND
                assert d.density == G.edge[u][v]['data'][id].density
                # subtract density from unstranded and add to stranded
                G.edge[u][v]['strand_density'] += d.density * redistribute_arr 

def create_strand_specific_graphs(G):
    '''
    build separate subgraphs of G for forward/reverse/unknown strands
   
    nodes in the graph will have a 'density' attribute equal to the
    total strand-specific density for that node

    edges (u,v) in the graph will have a 'density' attribute corresponding
    to the fraction of total coverage flowing out of node 'u' that goes 
    to 'v'
    
    NOTE: the unknown strand graph will not have any edges, just nodes
    
    returns a tuple containing (DiGraph,DiGraph,Graph) objects corresponding
    to the forward/reverse/unknown strand
    '''
    # make graphs
    GG = (nx.DiGraph(), nx.DiGraph(), nx.Graph())
    # add nodes
    for n,d in G.nodes_iter(data=True):
        strand_density = d['strand_density']
        for strand in (POS_STRAND, NEG_STRAND, NO_STRAND):
            density = strand_density[strand]
            # use a slightly non-zero value to overcome any numerical precision issues
            if density > 1e-8:
                GG[strand].add_node(n, density=density) 
    # add edges
    for u,v,d in G.edges_iter(data=True):            
        # determine strand based on coordinates of nodes
        strand = int(u.start >= v.end)        
        # check against strand attribute
        # TODO: can remove this
        strands = set(x.strand for x in d['data'].itervalues())
        assert all(cmp_strand(x, strand) for x in strands)
        # consider both the edge strand as well as unknown (non-stranded) data
        # and since non-stranded data uses a non-directional graph can add edges
        # in either direction and it is the same        
        for s in (strand, NO_STRAND):
            # get strand density at this edge
            density = d['strand_density'][s]
            # add edge if there is positive density
            # use a slightly non-zero value to overcome any numerical 
            # precision issues
            if density > 1e-8:
                GG[s].add_edge(u, v, density=density)
    return GG

def get_chains(G, criteria_func):
    """
    find chains of nodes that are genomically adjacent
    and have no other edges 
    """
    imin2 = lambda x,y: x if x<=y else y 
    imax2 = lambda x,y: x if x>=y else y 
    node_chain_map = {}
    chains = {}
    # initialize each node to be in a "chain" by itself
    for n in G.nodes_iter():
        node_chain_map[n] = n
        chains[n] = set((n,))
    for u,v in G.edges_iter():
        if not criteria_func(G,u,v):
            continue
        # get chains containing these nodes
        u_new = node_chain_map[u]
        u_chain = chains[u_new]
        del chains[u_new]
        v_new = node_chain_map[v]
        v_chain = chains[v_new]
        del chains[v_new]
        # merge chains        
        merged_chain = u_chain.union(v_chain)
        merged_node = Exon(imin2(u_new.start, v_new.start),
                           imax2(u_new.end, v_new.end))
        # point all nodes in chain to new parent
        for n in merged_chain:
            node_chain_map[n] = merged_node
        chains[merged_node] = merged_chain
    return node_chain_map, chains

def add_chain_nodes(H, G, chains):
    for parent, nodes in chains.iteritems():
        # calculate density of all nodes in chain
        total_mass = 0.0
        total_length = 0
        for n in nodes:
            length = (n.end - n.start)
            total_mass += length * G.node[n]['density']
            total_length += length
        density = total_mass / float(total_length)
        # sort nodes by genome position and find the min/max
        sorted_nodes = sorted(nodes, key=operator.attrgetter('start'))
        # make new node
        H.add_node(parent, chain=sorted_nodes, density=density, 
                   length=total_length)

def add_chain_edges(H, G, node_chain_map):
    for u,v,d in G.edges_iter(data=True):
        u_chain_node = node_chain_map[u]
        v_chain_node = node_chain_map[v]
        if u_chain_node != v_chain_node:
            H.add_edge(u_chain_node, v_chain_node, attr_dict=d)

def collapse_chains(G, directed=True):
    """
    find groups of nodes are genomically adjacent e.g. (100,200)-(200,300)
    and group them into a single node

    returns new DiGraph object.  each node has a 'chain' attribute 
    containing the child nodes making up the chain.  nodes also have
    'density' and 'length' attributes    
    """
    def directed_func(G,u,v):
        # check start/end coordinates
        if (u.start != v.end) and (u.end != v.start):
            return False
        # see if edge nodes have degree larger than '1'
        if ((G.out_degree(u) > 1) or (G.in_degree(v) > 1)):
            return False
        return True
    def undirected_func(G,u,v):
        # check start/end coordinates
        if (u.start != v.end) and (u.end != v.start):
            return False
        return True
    H = nx.DiGraph()
    if len(G) > 0:
        func = directed_func if directed else undirected_func
        node_chain_map, chains = get_chains(G, criteria_func=func)
        add_chain_nodes(H, G, chains)
        add_chain_edges(H, G, node_chain_map)
    return H

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
    G, transcript_graph_map, unstranded_ids = \
        add_transcripts_to_graph(transcripts, overhang_threshold)
    # compute fraction of fwd/rev strand coverage at each node
    # this creates a 'strand_density' array at each node/edge
    sum_strand_densities(G)
    # reallocate unstranded transcripts to fwd/rev strand 
    redistribute_unstranded_density(G, unstranded_ids, transcript_graph_map)
    # separate graph into strand-specific subgraphs
    HH = create_strand_specific_graphs(G)
    # collapse chains of nodes in subgraphs, which reformats the graph
    # attributes
    GG = (collapse_chains(HH[POS_STRAND]),
          collapse_chains(HH[NEG_STRAND]),
          collapse_chains(HH[NO_STRAND], directed=False))
    return GG
