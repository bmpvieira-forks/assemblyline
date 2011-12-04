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
from assembler import assemble_transcript_graph

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

def add_exon(exon, tdata, boundaries):
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
        exon_tdata = TranscriptData(id=tdata.id, strand=tdata.strand, 
                                    density=tdata.density)
        nodes.append((n, exon_tdata))
        # make edges between split exon according to 
        # strand being assembled.
        if len(nodes) > 0:
            if cmp_strand(tdata.strand, NEG_STRAND):
                edges.append((n, nodes[-1], tdata))
            if cmp_strand(tdata.strand, POS_STRAND):
                edges.append((nodes[-1], n, tdata))
    assert nodes[-1][0].end == exon.end
    return nodes, edges

def add_transcript(transcript, boundaries):
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
    enodes, eedges = add_exon(exons[0], tdata, boundaries)
    tnodes.extend(enodes)
    tedges.extend(eedges)
    for e2 in exons[1:]:
        # add exon
        enodes, eedges = add_exon(e2, tdata, boundaries)
        # add intron -> exon edges
        if strand != NO_STRAND:
            if strand == NEG_STRAND:
                tedges.append((enodes[0], tnodes[-1][0], tdata))
            else:
                tedges.append((tnodes[-1][0], enodes[0], tdata))
        # update nodes/edges
        tnodes.extend(enodes)
        tedges.extend(eedges)
    return tnodes, tedges

def sum_strand_densities(G):
    '''
    compute the fraction of (+) versus (-) transcription at each node
    '''
    nodes = G.nodes()
    for n in nodes:
        d = G.node[n]
        # sum density at each node
        strand_density = [0.0, 0.0, 0.0]
        for ndata in d['data']:
            strand_density[ndata.strand] += ndata.density
        d['strand_density'] = strand_density

def redistribute_unstranded_density(G, transcripts, transcript_nodes):
    '''
    build separate subgraphs of G - a forward strand graph and a 
    reverse strand graph (if necessary)
    
    pos_frac: fraction of all 'stranded' coverage that is on the 
    positive strand.  this is used to partition unstranded coverage
    proportionally
 
    nodes in the graph will have a 'weight' attribute equal to the
    total strand-specific coverage for that node

    edges (u,v) in the graph will have a 'weight' attribute corresponding
    to the fraction of total coverage flowing out of node 'u' that goes 
    to 'v'
    '''
    # compute fraction of fwd/rev strand coverage at each node    
    sum_strand_densities(G)
    # compute the fraction of fwd/rev coverage for each transcript
    for t in transcripts:
        # only re-apportion the density for unstranded transcripts
        if t.strand != NO_STRAND:
            continue
        tnodes = transcript_nodes[t.id]        
        # sum the coverage mass across the transcript
        strand_mass = np.zeros(3,float)
        for n in tnodes:
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
            strand_fracs = np.array((pos_frac, 1.0 - pos_frac, -1.0))
            for n in tnodes:
                # subtract density from unstranded and add to stranded
                G.node[n]['strand_density'] += (t.density * strand_fracs)




class TranscriptGraph(object):
    def __init__(self):
        pass

    def get_exon_ids(self, n):
        exon_data_list = self.G.node[n]['data']
        return [x.id for x in exon_data_list]

    def _add_node(self, n, tdata):
        if n not in self.G:  
            self.G.add_node(n, data=[])        
        nd = self.G.node[n]
        nd['data'].append(tdata)

    def _add_edge(self, u, v, tdata):
        if not self.G.has_edge(u, v):
            self.G.add_edge(u, v, data=[])
        ed = self.G.edge[u][v]
        ed['data'].append(tdata)
    
    def _add_exon(self, exon, tdata, boundaries):
        """
        exon: Exon object to split
        tdata: TranscriptData object with 'density' equal to the coverage 
        density
        """
        nodes = []
        edges = []
        for start, end in split_exon(exon, boundaries):
            n = Exon(start, end)
            exon_tdata = TranscriptData(id=tdata.id, strand=tdata.strand, 
                                        density=tdata.density)
            self._add_node(n, exon_tdata)
            # add edges between split exon according to 
            # strand being assembled.
            if len(nodes) > 0:
                if cmp_strand(tdata.strand, NEG_STRAND):
                    self._add_edge(n, nodes[-1], tdata)
                if cmp_strand(tdata.strand, POS_STRAND):
                    self._add_edge(nodes[-1], n, tdata)
            # continue loop
            nodes.append(n)
        assert nodes[-1].end == exon.end
        return nodes

    def _add_transcript(self, transcript, boundaries):
        exons = transcript.exons
        strand = transcript.strand
        tdata = TranscriptData(id=transcript.id, strand=strand, 
                               density=transcript.fpkm)
        # split exons that cross boundaries and get the
        # new exons of the path
        split_exons = []
        # add the first exon to initialize the loop
        # (all transcripts must have at least one exon)
        split_exons.extend(self._add_exon(exons[0], tdata, boundaries))
        for e2 in exons[1:]:
            # add exon
            e2_exons = self._add_exon(e2, tdata, boundaries)
            # add intron -> exon edges
            if strand != NO_STRAND:
                if strand == NEG_STRAND:
                    self._add_edge(e2_exons[0], split_exons[-1], tdata)
                else:
                    self._add_edge(split_exons[-1], e2_exons[0], tdata)
            # continue loop
            split_exons.extend(e2_exons)
        return split_exons

    @staticmethod
    def from_transcripts(transcripts, overhang_threshold=0):
        '''
        overhang_threshold: integer greater than zero specifying the 
        maximum exon overhang that can be trimmed to match an intron
        boundary.  exons that overhang more than this will be 
        considered independent transcript start sites or end sites
        
        note: this method cannot be called multiple times.  each time this
        function is invoked, the previously stored transcripts will be 
        deleted and overwritten        
        '''
        tg = TranscriptGraph()
        tg.G = nx.DiGraph()
        # trim the transcripts (modifies transcripts in place)
        trim_transcripts(transcripts, overhang_threshold)
        # find the intron domains of the transcripts
        boundaries = find_exon_boundaries(transcripts)
        # keep mapping of transcripts to graph nodes
        unstranded_transcripts = []
        # add transcripts
        for t in transcripts:
            # save NO_STRAND transcripts for later
            if t.strand == NO_STRAND:
                unstranded_transcripts.append(t)
            else:
                nodes = tg._add_transcript(t, boundaries)
        # now add unstranded transcripts
        for t in unstranded_transcripts:            
            nodes = tg._add_transcript(t, boundaries)

        # reallocate unstranded transcript coverage to fwd/rev strand 
        redistribute_unstranded_density(tg.G, transcripts, transcript_nodes)
        return tg

    def assemble(self, max_paths, fraction_major_path=0.10):
        if fraction_major_path <= 0:
            fraction_major_path = 1e-8        
        for res in assemble_transcript_graph(self.G, fraction_major_path, 
                                             max_paths):
            yield res