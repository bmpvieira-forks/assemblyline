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

from bx.intersection import Interval, IntervalTree
from bx.cluster import ClusterTree
from cnode import Node, cmp_strand, strand_int_to_str
from base import EXON, POS_STRAND, NEG_STRAND, NO_STRAND

class ExonData(object):
    __slots__ = ('id', 'strand', 'score')
    def __init__(self, id, strand, score):
        self.id = id
        self.strand = strand
        self.score = score
    def __repr__(self):
        return ("<%s(id='%s', strand='%s', score=%f')>" % 
                (self.__class__.__name__, self.id, self.strand, 
                 self.score)) 
    def __str__(self):
        return ("[ID=%s,cov=%.2f,%s]" % 
                (self.id, self.score, strand_int_to_str(self.strand)))

def find_transcript_boundaries(transcripts):
    '''
    input: a list of transcripts (not Node objects, these are transcripts)
    parsed directly from an input file and not added to an isoform graph
    yet. 
    
    output: sorted list of intron boundaries
    
    WARNING: multi-exonic unstranded transcripts will be ignored
    '''
    intron_boundaries = set()
    # first add introns to the graph and keep track of
    # all intron boundaries
    for transcript in transcripts:
        # add transcript exon boundaries
        for exon in transcript.exons:
            # keep track of positions where introns can be joined to exons
            intron_boundaries.add(exon.start)
            intron_boundaries.add(exon.end)
    # sort the intron boundary positions and add them to interval trees
    return sorted(intron_boundaries)

def get_in_degree(G, n, strand):
    d = 0
    for p in G.predecessors(n):
        print 'node', n, 'pred', p, 'strand', strand
        # ignore edges that are not compatible with the
        # strand being considered
        if (strand == POS_STRAND) and (p.end > n.start):
            continue
        elif (strand == NEG_STRAND) and (p.start < n.end):
            continue
        print 'YES'
        # found a valid predecessor
        d += 1
    return d

def get_out_degree(G, n, strand):
    d = 0
    for s in G.successors(n):
        # ignore edges that are not compatible with the
        # strand being considered
        if (strand == POS_STRAND) and (s.start < n.end):
            continue
        elif (strand == NEG_STRAND) and (s.end > n.start):
            continue
        # found a valid predecessor
        d += 1
    return d

def find_start_and_end_nodes(G, strand):
    # find unique starting positions and their 
    # corresponding nodes in the graph
    tss_node_dict = collections.defaultdict(lambda: [])
    end_nodes = []
    for n in G.nodes_iter():
        # look for start nodes
        in_degree = get_in_degree(G, n, strand)
        if in_degree == 0:
            tss_pos = n.end if strand == NEG_STRAND else n.start
            tss_node_dict[tss_pos].append(n)
        # look for end nodes
        out_degree = get_out_degree(G, n, strand)
        if out_degree == 0:
            end_nodes.append(n)
    logging.debug("START NODES (TSSs): %s" % tss_node_dict.values())
    logging.debug("END NODES: %s" % end_nodes)
    return tss_node_dict.values(), end_nodes


def assemble(G, strand):
    # add dummy start and end nodes to the graph
    start_node_lists, end_nodes = find_start_and_end_nodes(G, strand)    
    H = G.copy()
    # separate plus and minus strand assemblies
    # make a copy of the graph

    # delete edges from the opposite strand
    for edge in H.edges_iter():

        print edge
    
#    # add exon-exon edges
#    for n,attr_dict in G.nodes_iter(data=True):
#        if n.node_type == INTRON:
#            weight = sum(id_score_map[id] for id in attr_dict['ids'])            
#            for pred in G.predecessors_iter(n):
#                for succ in G.successors_iter(n):   
#                    psi = (weight/H.node[pred]['weight'])               
#                    H.add_edge(pred, succ, weight=weight, 
#                               psi=psi, ids=attr_dict['ids'])
#    return H
#    for Hsubgraph in nx.weakly_connected_component_subgraphs(H):
#        # add dummy start and end nodes to the graph
#        start_node_lists, end_nodes = find_start_and_end_nodes(Hsubgraph)    


class IsoformGraph(object):
    def __init__(self):
        self.G = nx.DiGraph()

    @staticmethod
    def from_transcripts(transcripts):
        g = IsoformGraph()
        g.add_transcripts(transcripts)

    def _add_exon_node(self, start, end, transcript):
        n = Node(start, end, NO_STRAND, EXON)
        if n not in self.G:  
            self.G.add_node(n, data=[])        
        nd = self.G.node[n]
        score = transcript.score * float(end - start) / transcript.length
        nd['data'].append(ExonData(id=transcript.id, strand=transcript.strand, score=score))
        return n

    def _add_exon(self, exon, transcript, boundaries):
        #print 'adding exon=%s' % str(exon)        
        strand = transcript.strand
        # find the indexes into the intron boundaries list that
        # border the exon.  all the indexes in between these two
        # are overlapping the exon and we must use them to break
        # the exon into pieces 
        start_ind = bisect.bisect_right(boundaries, exon.start)
        end_ind = bisect.bisect_left(boundaries, exon.end)
        exon_splits = [exon.start] + boundaries[start_ind:end_ind] + [exon.end]
        first_node = None
        prev_node = None
        for j in xrange(1, len(exon_splits)):            
            cur_node = self._add_exon_node(exon_splits[j-1], exon_splits[j], transcript)
            # add edges between split exon according to 
            # strand being assembled.  this allows edges between
            # split single exons (unstranded) to be oriented 
            # correctly
            if prev_node is None:
                first_node = cur_node
            else:
                if cmp_strand(strand, NEG_STRAND):
                    self.G.add_edge(cur_node, prev_node)
                if cmp_strand(strand, POS_STRAND):
                    self.G.add_edge(prev_node, cur_node)
            # continue loop
            prev_node = cur_node
        assert cur_node.end == exon.end
        return first_node, cur_node

    def _add_transcript(self, transcript, boundaries):
        #print 'adding transcript=%s' % transcript
        exons = transcript.exons        
        # add the first exon to initialize the loop
        # (all transcripts must have at least one exon)
        e1_start_node, e1_end_node = \
            self._add_exon(exons[0], transcript, boundaries)
        for e2 in exons[1:]:
            # add exon
            e2_start_node, e2_end_node = \
                self._add_exon(e2, transcript, boundaries)
            # unstranded transcripts may have multiple pieces
            # but will not create 'intron' edges between them
            if transcript.strand != NO_STRAND:            
                # add edges from exon -> exon
                if transcript.strand == NEG_STRAND:
                    self.G.add_edge(e2_start_node, e1_end_node)
                else:
                    self.G.add_edge(e1_end_node, e2_start_node)
            # continue loop
            e1_end_node = e2_end_node

    def add_transcripts(self, transcripts):
        # find the intron domains of the transcripts
        boundaries = find_transcript_boundaries(transcripts)
        for t in transcripts:
            self._add_transcript(t, boundaries)

    def trim(self):
        pass

    def assemble(self):
        # get connected components        
        for Gsub in nx.weakly_connected_component_subgraphs(self.G):
            assemble(Gsub, POS_STRAND)
