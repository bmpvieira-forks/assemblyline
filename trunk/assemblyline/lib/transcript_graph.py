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
    __slots__ = ('id', 'strand', 'scores')
    def __init__(self, id, strand, scores):
        self.id = id
        self.strand = strand
        self.scores = map(float, scores)
    def __repr__(self):
        return ("<%s(id='%s',strand='%s',scores='%s')>" % 
                (self.__class__.__name__, self.id, self.strand, 
                 str(self.scores))) 
    def __str__(self):
        return ("[id=%s,strand=%s,scores=%s]" % 
                (self.id, strand_int_to_str(self.strand),
                 '(%.2f,%.2f)' % (self.scores[POS_STRAND], 
                                  self.scores[NEG_STRAND])))

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

def calculate_strand_fraction(G, exon_splits):
    strand_scores = [0, 0]
    for start, end, score in exon_splits:
        n = Exon(start, end)
        if n not in G:
            continue
        exon_data_list = G.node[n]['data']
        for edata in exon_data_list:
            assert edata.strand != NO_STRAND
            if edata.strand == POS_STRAND:
                assert edata.scores[POS_STRAND] > 0
                assert edata.scores[NEG_STRAND] == 0
            elif edata.strand == NEG_STRAND:
                assert edata.scores[POS_STRAND] == 0
                assert edata.scores[NEG_STRAND] > 0                
            strand_scores[POS_STRAND] += edata.scores[POS_STRAND]
            strand_scores[NEG_STRAND] += edata.scores[NEG_STRAND]
    total_score = strand_scores[POS_STRAND] + strand_scores[NEG_STRAND]
    # if there is no "stranded" coverage at this node, then assign coverage
    # the positive strand by convention
    if total_score == 0:
        return 1.0
    # proportionally assign unstranded coverage based on amount of
    # plus and minus strand coverage
    pos_frac = strand_scores[POS_STRAND] / float(total_score)
    return pos_frac

def add_unstranded_transcript(G, transcript, boundaries, nbunch, ebunch):
    strand = transcript.strand
    cov = transcript.score / transcript.length
    for exon in transcript.exons:
        exon_splits = list(split_exon(exon, cov, boundaries))
        # the coverage allocation is computed here
        pos_frac = calculate_strand_fraction(G, exon_splits)
        prev_node = None
        # create nodes and edges to add to the graph
        for start, end, score in exon_splits:             
            scores = pos_frac * score, (1.0 - pos_frac) * score                    
            cur_node = Exon(start, end)
            exon_data = ExonData(id=transcript.id, strand=strand, scores=scores)
            nbunch.append((cur_node, exon_data))                    
            if prev_node is not None:
                # only add edges if there is evidence of coverage on the 
                # appropriate strand
                if pos_frac > 0:
                    ebunch.append((prev_node, cur_node, transcript.id, strand, scores[POS_STRAND]))
                if pos_frac < 1:
                    ebunch.append((cur_node, prev_node, transcript.id, strand, scores[NEG_STRAND]))
            # continue loop
            prev_node = cur_node

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
        ed['data'].append(EdgeData(id=id, score=score))
    
    def _add_stranded_exon(self, exon, id, strand, cov, boundaries):
        assert strand != NO_STRAND        
        nfirst, ndatafirst = None, None
        n1, ndata1 = None, None
        n2, ndata2 = None, None
        for start, end, score in split_exon(exon, cov, boundaries):
            if strand == POS_STRAND:
                scores = (score, 0)
            elif strand == NEG_STRAND:
                scores = (0, score)
            n2 = Exon(start, end)
            ndata2 = ExonData(id=id, strand=strand, scores=scores)
            self._add_node(n2, ndata2)
            # add edges between split exon according to 
            # strand being assembled.  this allows edges between
            # split single exons (unstranded) to be oriented 
            # correctly
            if n1 is None:
                nfirst = n2
                ndatafirst = ndata2
            elif strand == NEG_STRAND:
                self._add_edge(n2, n1, id, strand, ndata2.scores[NEG_STRAND])
            else:
                assert strand == POS_STRAND
                self._add_edge(n1, n2, id, strand, ndata1.scores[POS_STRAND])
            # continue loop
            n1, ndata1 = n2, ndata2
        assert n2.end == exon.end
        return nfirst, ndatafirst, n2, ndata2

    def _add_stranded_transcript(self, transcript, boundaries):
        assert transcript.strand != NO_STRAND
        exons = transcript.exons        
        strand = transcript.strand
        cov = transcript.score / transcript.length
        # add the first exon to initialize the loop
        # (all transcripts must have at least one exon)
        e1_start_node, e1_start_data, e1_end_node, e1_end_data = \
            self._add_stranded_exon(exons[0], transcript.id, strand, cov, boundaries)
        for e2 in exons[1:]:
            # add exon
            e2_start_node, e2_start_data, e2_end_node, e2_end_data = \
                self._add_stranded_exon(e2, transcript.id, strand, cov, boundaries)
            # add edges from exon -> exon
            if strand == NEG_STRAND:
                self._add_edge(e2_start_node, e1_end_node, transcript.id, 
                               strand, e2_start_data.scores[NEG_STRAND])
            else:
                self._add_edge(e1_end_node, e2_start_node, transcript.id,
                               strand, e1_end_data.scores[POS_STRAND])
            # continue loop
            e1_end_node = e2_end_node
            e1_end_data = e2_end_data

    def _add_unstranded_transcripts(self, transcripts, boundaries):
        nbunch = []
        ebunch = []
        for transcript in transcripts:
            add_unstranded_transcript(self.G, transcript, boundaries, nbunch, ebunch)
        # now add all the nodes and edges
        for n,exon_data in nbunch:
            self._add_node(n, exon_data)
        for u,v,id,strand,score in ebunch:
            self._add_edge(u, v, id, strand, score)

    def add_transcripts(self, transcripts):
        '''
        note: this method cannot be called multiple times.  each time this
        function is invoked, the previously stored transcripts will be 
        deleting and overwritten
        '''
        self.G = nx.DiGraph()
        # find the intron domains of the transcripts
        boundaries = find_transcript_boundaries(transcripts)
        # add the stranded transcripts first and save the
        # unstranded transcripts for later
        nostrand_transcripts = []
        for t in transcripts:
            if t.strand == NO_STRAND:
                nostrand_transcripts.append(t)
                continue
            self._add_stranded_transcript(t, boundaries)
        # add the unstranded transcripts
        self._add_unstranded_transcripts(nostrand_transcripts, boundaries)

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
