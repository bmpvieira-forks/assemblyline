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
from base import Exon, POS_STRAND, NEG_STRAND, NO_STRAND, cmp_strand, strand_int_to_str

from assembler import assemble_transcript_graph

class ExonData(object):
    __slots__ = ('id', 'strand', 'score')
    def __init__(self, id, strand, score):
        self.id = id
        self.strand = strand
        self.score = float(score)
    def __repr__(self):
        return ("<%s(id='%s',strand='%s',score='%s')>" % 
                (self.__class__.__name__, self.id, self.strand, 
                 str(self.score))) 
    def __str__(self):
        return ("[id=%s,strand=%s,score=%.2f]" % 
                (self.id, strand_int_to_str(self.strand), self.score))

class EdgeData(object):
    __slots__ = ('id', 'score')
    def __init__(self, id, score):
        self.id = id
        self.score = float(score)
    def __repr__(self):
        return ("<%s(id='%s', score='%f')>" % 
                (self.__class__.__name__, self.id, self.score)) 
    def __str__(self):
        return ("[id=%s,score=%.2f]" % (self.id, self.score))

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

def get_free_ends(strand, in_degree, out_degree):
    """returns whether the (left,right) genomic ends of the node are free ends without edges"""
    # TODO: move to C code
    fiveprime = (in_degree == 0)
    threeprime = (out_degree == 0)
    if strand == NEG_STRAND:
        return threeprime, fiveprime
    else:
        return fiveprime, threeprime

def trim_interval(start, end, trim_left, trim_right, 
                  intron_starts, intron_ends,
                  overhang_threshold):
    trim_start, trim_end = start, end
    if trim_left:
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
    # if both start and end were trimmed it is
    # possible that they could be trimmed to 
    # match different introns and generate a 
    # weird and useless trimmed exon that is not
    # compatible with either the left or right
    # introns. resolve these situations by choosing
    # the smaller of the two trimmed distances as
    # the more likely
    if trim_end < trim_start:
        left_trim_dist = trim_start - start
        right_trim_dist = end - trim_end
        if left_trim_dist <= right_trim_dist:
            trim_start, trim_end = trim_start, end
        else:
            trim_start, trim_end = start, trim_end
    return trim_start, trim_end


def split_exon(exon, transcript_coverage, boundaries):
    # find the indexes into the intron boundaries list that
    # border the exon.  all the indexes in between these two
    # are overlapping the exon and we must use them to break
    # the exon into pieces 
    start_ind = bisect.bisect_right(boundaries, exon.start)
    end_ind = bisect.bisect_left(boundaries, exon.end)
    exon_splits = [exon.start] + boundaries[start_ind:end_ind] + [exon.end]
    for j in xrange(1, len(exon_splits)):
        start, end = exon_splits[j-1], exon_splits[j]
        score = float(end - start) * transcript_coverage
        yield start, end, score
        
class TranscriptGraph(object):
    def __init__(self):
        pass

    @staticmethod
    def from_transcripts(transcripts):
        g = TranscriptGraph()
        g.add_transcripts(transcripts)
        return g

    def _add_node(self, n, ndata):
        if n not in self.G:  
            self.G.add_node(n, data=[])        
        nd = self.G.node[n]
        nd['data'].append(ndata)

    def _add_edge(self, u, v, id, strand, score):        
        if not self.G.has_edge(u, v):
            self.G.add_edge(u, v, strand=strand, data=[])
        ed = self.G.edge[u][v]
        assert ed['strand'] == strand 
        ed['data'].append(EdgeData(id=id, score=score))
    
    def _add_exon(self, exon, id, strand, cov, boundaries):
        assert strand != NO_STRAND        
        nfirst, ndatafirst = None, None
        n1, ndata1 = None, None
        n2, ndata2 = None, None
        for start, end, score in split_exon(exon, cov, boundaries):
            n2 = Exon(start, end)
            ndata2 = ExonData(id=id, strand=strand, score=score)
            self._add_node(n2, ndata2)
            # add edges between split exon according to 
            # strand being assembled.
            if n1 is None:
                nfirst = n2
                ndatafirst = ndata2
            else:
                if cmp_strand(strand, NEG_STRAND):
                    self._add_edge(n2, n1, id, strand, score)
                if cmp_strand(strand, POS_STRAND):
                    self._add_edge(n1, n2, id, strand, score)
            # continue loop
            n1, ndata1 = n2, ndata2
        assert n2.end == exon.end
        return nfirst, ndatafirst.score, n2, ndata2.score

    def _add_transcript(self, transcript, boundaries):
        assert transcript.strand != NO_STRAND
        exons = transcript.exons        
        strand = transcript.strand
        cov = transcript.score / transcript.length
        # add the first exon to initialize the loop
        # (all transcripts must have at least one exon)
        e1_start_node, e1_start_score, e1_end_node, e1_end_score = \
            self._add_exon(exons[0], transcript.id, strand, cov, boundaries)
        for e2 in exons[1:]:
            # add exon
            e2_start_node, e2_start_score, e2_end_node, e2_end_score = \
                self._add_exon(e2, transcript.id, strand, cov, boundaries)
            # add intron -> exon edges
            if strand != NO_STRAND:
                if strand == NEG_STRAND:
                    self._add_edge(e2_start_node, e1_end_node, transcript.id, 
                                   strand, e2_start_score)
                else:
                    self._add_edge(e1_end_node, e2_start_node, transcript.id,
                                   strand, e1_end_score)
            # continue loop
            e1_end_node = e2_end_node
            e1_end_score = e2_end_score

    def add_transcripts(self, transcripts):
        '''
        note: this method cannot be called multiple times.  each time this
        function is invoked, the previously stored transcripts will be 
        deleting and overwritten
        '''
        self.G = nx.DiGraph()
        # find the intron domains of the transcripts
        boundaries = find_transcript_boundaries(transcripts)
        # add transcripts
        for t in transcripts:
            self._add_transcript(t, boundaries)

    def get_exon_ids(self, n):
        exon_data_list = self.G.node[n]['data']
        return [x.id for x in exon_data_list]

    def trim_introns(self):
        pass
    
    def trim_exons(self):
        # trim exons before computing exon domains
        if overhang_threshold > 0:
            left_free, right_free = get_free_ends(exon.strand, len(preds), len(succs))
            start, end = trim_interval(exon.start, exon.end, 
                                       left_free, right_free,
                                       intron_starts, intron_ends,
                                       overhang_threshold)        

    def assemble(self, max_paths, fraction_major_isoform=0.10):
        for res in assemble_transcript_graph(self.G, fraction_major_isoform, 
                                             max_paths):
            yield res
