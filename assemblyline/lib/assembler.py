'''
Created on Feb 20, 2011

@author: mkiyer
'''
import collections
import logging
import networkx as nx

from bx.cluster import ClusterTree

from base import Exon, POS_STRAND, NEG_STRAND, NO_STRAND, cmp_strand

# constants for fake 'start' and 'end' nodes
DUMMY_START_NODE = Exon(-1,-1)
DUMMY_END_NODE = Exon(-2,-2)

# constants for use as graph attributes
NODE_WEIGHT = 'weight'
NODE_TSS_ID = 'tss_id'

EDGE_WEIGHT = 'weight'
EDGE_DENSITY = 'density'
EDGE_OUT_FRAC = 'outfrac'
EDGE_IN_FRAC = 'infrac'

PLEN = 'plen'
PWEIGHT = 'pweight'
PDENSITY = 'pdensity'
PSRC = 'psrc'


#def calculate_strand_fraction(G):
#
#    strand_scores = [0, 0, 0]
#    for n,d in G.nodes_iter(data=True):
#        for ndata in d['data']:
#            strand_scores[ndata.strand] += ndata.score / float((n.end - n.start))        
#    total_score = strand_scores[POS_STRAND] + strand_scores[NEG_STRAND]
#    # if there is no "stranded" coverage at this node, then assign coverage
#    # the positive strand by convention
#    if total_score == 0:
#        return 1.0
#    # proportionally assign unstranded coverage based on amount of
#    # plus and minus strand coverage
#    pos_frac = strand_scores[POS_STRAND] / float(total_score)
#    return pos_frac

def calculate_strand_fraction(G):
    '''
    compute the fraction of (+) versus (-) transcription in the
    graph by using the summed coverage densities in all the nodes 
    '''
    cluster_tree = ClusterTree(0,1)
    nodes = G.nodes()
    for i,n in enumerate(nodes):    
        cluster_tree.insert(n.start, n.end, i)
    for start, end, indexes in cluster_tree.getregions():
        strand_scores = [0, 0, 0]
        for i in indexes:
            n = nodes[i]
            for ndata in G.node[n]['data']:
                strand_scores[ndata.strand] += ndata.score / float((n.end - n.start))        
        total_score = strand_scores[POS_STRAND] + strand_scores[NEG_STRAND]        
        # if there is no "stranded" coverage at this node, then assign coverage
        # the positive strand by convention
        if total_score == 0:
            return 1.0
        # proportionally assign unstranded coverage based on amount of
        # plus and minus strand coverage
        pos_frac = strand_scores[POS_STRAND] / float(total_score)
        yield G.subgraph(nodes), pos_frac
    del cluster_tree

def sum_node_weights(exon_data_list):
    '''
    sum exon data scores by strand
    '''
    # sum coverage of node by strand
    strand_weights = [0, 0, 0]
    for edata in exon_data_list:
        strand_weights[edata.strand] += edata.score
    return strand_weights

def build_strand_specific_graphs(G, pos_frac):
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
    # make graphs
    GG = (nx.DiGraph(), nx.DiGraph())
    # cluster nodes and compute fraction of fwd/rev strand coverage    
    for Gsub, pos_frac in calculate_strand_fraction(G.nodes()):
        # add nodes first
        strand_fracs = (pos_frac, 1.0 - pos_frac)
        for n in Gsub.nodes_iter():
            strand_weights = sum_node_weights(G.node[n]['data'])
            # partition the unstranded coverage to fwd/rev strands
            # according to the fraction of stranded coverage observed
            for strand in (POS_STRAND, NEG_STRAND):
                w = strand_weights[strand] + strand_fracs[strand] * strand_weights[NO_STRAND] 
                if w > 0:
                    GG[strand].add_node(n, weight=w)
        # add edges
        for u,v,d in G.edges_iter(data=True):
            # determine strand based on coordinates of nodes
            strand = int(u.start >= v.end)
            # check against strand attribute
            strandattr = d['strand']
            assert cmp_strand(strand, strandattr)
            # get total coverage on this edge
            w = sum(x.score for x in d['data'])
            if strandattr == NO_STRAND:
                w *= strand_fracs[strand]
            if w > 0:
                GG[strand].add_edge(u, v, attr_dict={EDGE_WEIGHT: w})   
    return GG


#def add_unstranded_transcript(G, transcript, boundaries, nbunch, ebunch):
#    strand = transcript.strand
#    cov = transcript.score / transcript.length
#    for exon in transcript.exons:
#        exon_splits = list(split_exon(exon, cov, boundaries))
#        # the coverage allocation is computed here
#        pos_frac = calculate_strand_fraction(G, exon_splits)
#        prev_node = None
#        # create nodes and edges to add to the graph
#        for start, end, score in exon_splits:             
#            scores = pos_frac * score, (1.0 - pos_frac) * score                    
#            cur_node = Exon(start, end)
#            exon_data = ExonData(id=transcript.id, strand=strand, scores=scores)
#            nbunch.append((cur_node, exon_data))                    
#            if prev_node is not None:
#                # only add edges if there is evidence of coverage on the 
#                # appropriate strand
#                if pos_frac > 0:
#                    ebunch.append((prev_node, cur_node, transcript.id, strand, scores[POS_STRAND]))
#                if pos_frac < 1:
#                    ebunch.append((cur_node, prev_node, transcript.id, strand, scores[NEG_STRAND]))
#            # continue loop
#            prev_node = cur_node
#    def _add_stranded_exon(self, exon, id, strand, cov, boundaries):
#        assert strand != NO_STRAND        
#        nfirst, ndatafirst = None, None
#        n1, ndata1 = None, None
#        n2, ndata2 = None, None
#        for start, end, score in split_exon(exon, cov, boundaries):
#            if strand == POS_STRAND:
#                scores = (score, 0)
#            elif strand == NEG_STRAND:
#                scores = (0, score)
#            n2 = Exon(start, end)
#            ndata2 = ExonData(id=id, strand=strand, scores=scores)
#            self._add_node(n2, ndata2)
#            # add edges between split exon according to 
#            # strand being assembled.  this allows edges between
#            # split single exons (unstranded) to be oriented 
#            # correctly
#            if n1 is None:
#                nfirst = n2
#                ndatafirst = ndata2
#            elif strand == NEG_STRAND:
#                self._add_edge(n2, n1, id, strand, ndata2.scores[NEG_STRAND])
#            else:
#                assert strand == POS_STRAND
#                self._add_edge(n1, n2, id, strand, ndata1.scores[POS_STRAND])
#            # continue loop
#            n1, ndata1 = n2, ndata2
#        assert n2.end == exon.end
#        return nfirst, ndatafirst, n2, ndata2
#
#    def _add_stranded_transcript(self, transcript, boundaries):
#        assert transcript.strand != NO_STRAND
#        exons = transcript.exons        
#        strand = transcript.strand
#        cov = transcript.score / transcript.length
#        # add the first exon to initialize the loop
#        # (all transcripts must have at least one exon)
#        e1_start_node, e1_start_data, e1_end_node, e1_end_data = \
#            self._add_stranded_exon(exons[0], transcript.id, strand, cov, boundaries)
#        for e2 in exons[1:]:
#            # add exon
#            e2_start_node, e2_start_data, e2_end_node, e2_end_data = \
#                self._add_stranded_exon(e2, transcript.id, strand, cov, boundaries)
#            # add edges from exon -> exon
#            if strand == NEG_STRAND:
#                self._add_edge(e2_start_node, e1_end_node, transcript.id, 
#                               strand, e2_start_data.scores[NEG_STRAND])
#            else:
#                self._add_edge(e1_end_node, e2_start_node, transcript.id,
#                               strand, e1_end_data.scores[POS_STRAND])
#            # continue loop
#            e1_end_node = e2_end_node
#            e1_end_data = e2_end_data
#    def _add_unstranded_transcripts(self, transcripts, boundaries):
#        nbunch = []
#        ebunch = []
#        for transcript in transcripts:
#            add_unstranded_transcript(self.G, transcript, boundaries, nbunch, ebunch)
#        # now add all the nodes and edges
#        for n,exon_data in nbunch:
#            self._add_node(n, exon_data)
#        for u,v,id,strand,score in ebunch:
#            self._add_edge(u, v, id, strand, score)
#
#    def add_transcripts(self, transcripts):
#        '''
#        note: this method cannot be called multiple times.  each time this
#        function is invoked, the previously stored transcripts will be 
#        deleting and overwritten
#        '''
#        self.G = nx.DiGraph()
#        # find the intron domains of the transcripts
#        boundaries = find_transcript_boundaries(transcripts)
#        # add the stranded transcripts first and save the
#        # unstranded transcripts for later
#        nostrand_transcripts = []
#        for t in transcripts:
#            if t.strand == NO_STRAND:
#                nostrand_transcripts.append(t)
#                continue
#            self._add_stranded_transcript(t, boundaries)
#        # add the unstranded transcripts
#        self._add_unstranded_transcripts(nostrand_transcripts, boundaries)
#            def _add_unstranded_transcripts(self, transcripts, boundaries):
#        nbunch = []
#        ebunch = []
#        for transcript in transcripts:
#            add_unstranded_transcript(self.G, transcript, boundaries, nbunch, ebunch)
#        # now add all the nodes and edges
#        for n,exon_data in nbunch:
#            self._add_node(n, exon_data)
#        for u,v,id,strand,score in ebunch:
#            self._add_edge(u, v, id, strand, score)



def find_start_and_end_nodes(G, strand):
    # find unique starting positions and their 
    # corresponding nodes in the graph
    tss_node_dict = collections.defaultdict(lambda: [])
    for n,d in G.in_degree_iter():
        if d == 0:
            tss_pos = n.end if strand == NEG_STRAND else n.start
            tss_node_dict[tss_pos].append(n)
    # annotate each TSS with a unique id
    tss_id = 0
    start_nodes = []
    for tss_nodes in tss_node_dict.itervalues():
        for n in tss_nodes:
            G.node[n][NODE_TSS_ID] = tss_id
            start_nodes.append(n)
        tss_id += 1
    end_nodes = [n for (n,d) in G.out_degree_iter()
                 if (d == 0)]
    return start_nodes, end_nodes

def visit_node_from_parent(G, parent, child):    
    # compute length, coverage, and rpkm of the current
    # path when extended to include the child
    pattrs = G.node[parent]
    # new path length is parent length + child length    
    path_length = pattrs.get(PLEN, 0) + (child.end - child.start)    
    # to compute the coverage going from the parent to the child, take
    # the existing path coverage at the parent node and multiply it by
    # the outgoing fraction, then add the child coverage multiplied by
    # its incoming fraction.
    parent_cov = pattrs.get(PWEIGHT, pattrs[NODE_WEIGHT])
    child_cov = G.node[child].get(NODE_WEIGHT, 0)
    edata = G.edge[parent][child]    
    path_cov = ((parent_cov * edata[EDGE_OUT_FRAC]) + 
                (child_cov * edata[EDGE_IN_FRAC]))
    # score is total coverage divided by the length of the path
    path_density = path_cov / path_length
    return path_length, path_cov, path_density

def dyn_prog_search(G, source):
    """Find the highest scoring path by dynamic programming"""
    # Taken from NetworkX source code http://networkx.lanl.gov
    stack = [source]
    while stack:
        parent = stack[0]
        for child in G.successors_iter(parent):
            logging.debug('Visiting %s -> %s' % (parent, child))            
            path_length, path_cov, path_density = visit_node_from_parent(G, parent, child)
            logging.debug('\tlen=%d cov=%f density=%f child_density=%f' %
                          (path_length, path_cov, path_density, G.node[child].get(PDENSITY,-1)))
            # only update this node if it has not been visited
            # or if the new score is higher than the old score
            cattrs = G.node[child]
            if path_density > cattrs.get(PDENSITY, -1):
                # keep pointer to parent node that produced this high scoring path
                cattrs[PLEN] = path_length
                cattrs[PWEIGHT] = path_cov
                cattrs[PDENSITY] = path_density
                cattrs[PSRC] = parent
                logging.debug('\tupdated child len=%d cov=%f density=%f' % 
                              (path_length, path_cov, path_density))
                # continue iterating through the child path
                stack.append(child)
        stack.pop(0)

def find_best_path(G, source, sink):
    '''
    use dynamic programming to find the highest scoring path through 
    the graph starting from 'source' and ending at 'sink'
    '''
    # run the dynamic programming algorithm to
    # score paths
    dyn_prog_search(G, source)
    # traceback
    path = [sink]
    pattrs = G.node[sink]
    score = pattrs[PDENSITY]
    length = pattrs[PLEN]
    cov = pattrs[PWEIGHT]
    while path[-1] != source:
        #logging.debug('path=%s node=%s attrs=%s' %
        #              (path, path[-1], G.node[path[-1]]))
        path.append(G.node[path[-1]][PSRC])
    path.reverse()
    #logging.debug("FINAL path=%s" % (path))
    return path, length, cov, score

def clear_path_attributes(G):
    '''
    remove dynamic programming node attributes that are 
    added to the graph temporarily
    
    must call this before calling the dynamic programming
    algorithm again on the same graph    
    '''
    # clear path attributes
    for n in G.nodes_iter():
        nattrs = G.node[n]
        if PLEN in nattrs:
            del nattrs[PLEN]
            del nattrs[PWEIGHT]
            del nattrs[PDENSITY]
            del nattrs[PSRC]    

def recalculate_edge_attrs(G, u, v):
    # find total weight leaving node 'u'
    u_succs = G.successors(u)
    out_weight_total = sum(G.edge[u][x][EDGE_WEIGHT] for x in u_succs)
    if out_weight_total == 0:
        out_weight_total = 1
    # compute edge attributes leaving node 'u'
    for succ in u_succs:
        eattrs = G.edge[u][succ]
        out_frac = eattrs[EDGE_WEIGHT]/float(out_weight_total)
        if out_frac == 0:
            density = 0
        else:
            density = (G.node[u][NODE_WEIGHT] * out_frac)/float(u.end - u.start)
        eattrs[EDGE_OUT_FRAC] = out_frac
        eattrs[EDGE_DENSITY] = density
    # find total density entering node 'v'
    v_preds = G.predecessors(v)
    in_density_total = sum(G.edge[x][v][EDGE_DENSITY] for x in v_preds)
    if in_density_total == 0:
        in_density_total = 1
    # compute edge density fractions entering node 'v'
    for pred in v_preds:
        eattrs = G.edge[pred][v]
        in_frac = eattrs[EDGE_DENSITY]/float(in_density_total)
        eattrs[EDGE_IN_FRAC] = in_frac

def subtract_path(G, path, path_density):
    '''
    subtract the coverage density incurred by traversing the 
    given path from each node and edge in the graph
    '''
    for i,v in enumerate(path):
        # subtract coverage from node
        v_cov = path_density * (v.end - v.start)
        G.node[v][NODE_WEIGHT] = max(0, G.node[v][NODE_WEIGHT] - v_cov)
        if i == 0:
            continue
        # subtract coverage from edge
        u = path[i-1]
        u_cov = path_density * (u.end - u.start)
        eattrs = G.edge[u][v]
        eattrs[EDGE_WEIGHT] = max(0, eattrs[EDGE_WEIGHT] - u_cov)
        # recompute edge attrs
        recalculate_edge_attrs(G, u, v)

def add_dummy_start_end_nodes(G, start_nodes, end_nodes):
    # add 'dummy' tss nodes if necessary
    G.add_node(DUMMY_START_NODE, weight=0)
    for start_node in start_nodes:        
        logging.debug('adding dummy %s -> %s' % (DUMMY_START_NODE, start_node))
        attr_dict = {EDGE_WEIGHT: 0.0,
                     EDGE_DENSITY: 0.0,
                     EDGE_OUT_FRAC: 0.0,
                     EDGE_IN_FRAC: 1.0}                     
        G.add_edge(DUMMY_START_NODE, start_node, attr_dict=attr_dict)
    # add a single 'dummy' end node
    G.add_node(DUMMY_END_NODE, weight=0)
    for end_node in end_nodes:
        logging.debug('adding dummy %s -> %s' % (end_node, DUMMY_END_NODE))
        attr_dict = {EDGE_WEIGHT: 0.0,
                     EDGE_DENSITY: 0.0,
                     EDGE_OUT_FRAC: 1.0,
                     EDGE_IN_FRAC: 0.0}
        G.add_edge(end_node, DUMMY_END_NODE, attr_dict=attr_dict)

def calculate_edge_attrs(G):
    for u,nbrdict in G.adjacency_iter():
        # find total coverage flowing out of this node
        out_score_total = sum(eattrs[EDGE_WEIGHT] for eattrs in nbrdict.itervalues())        
        for v, eattrs in nbrdict.iteritems():
            # find fraction flowing out of each node
            out_frac = eattrs[EDGE_WEIGHT]/float(out_score_total)
            # compute coverage density flowing out of edge
            weight = (G.node[u][NODE_WEIGHT] * out_frac)/float(u.end - u.start)
            eattrs[EDGE_OUT_FRAC] = out_frac
            eattrs[EDGE_DENSITY] = weight
    # now fraction of weight flowing in to each node
    # reverse edge directions and use adjacency iter again
    G.reverse(copy=False)
    for u,nbrdict in G.adjacency_iter():
        # find total coverage density flowing into this node
        in_weight_total = sum(eattrs[EDGE_DENSITY] for eattrs in nbrdict.itervalues())
        for v, eattrs in nbrdict.iteritems():
            # find fraction flowing into node            
            in_frac = eattrs[EDGE_DENSITY]/float(in_weight_total)
            eattrs[EDGE_IN_FRAC] = in_frac
    # reverse the edges back to normal
    G.reverse(copy=False)

def find_suboptimal_paths(G, start_nodes, end_nodes, 
                          fraction_major_path,
                          max_paths):
    # compute initial edge weights and coverage flow in/out of nodes
    calculate_edge_attrs(G)    
    # at an artificial node at the start and end so that all start
    # nodes are searched together for the best path
    add_dummy_start_end_nodes(G, start_nodes, end_nodes)
    score_limit = -1
    while True:        
        # find highest scoring path through graph
        path, path_length, path_cov, path_density = \
            find_best_path(G, DUMMY_START_NODE, DUMMY_END_NODE)
        logging.debug("Path density=%f limit=%f cov=%.2f length=%d path=%s" % 
                      (path_density, score_limit, path_cov, path_length, path))
        # set score limit if this is the first path
        if score_limit < 0:
            score_limit = path_density * fraction_major_path
        elif path_density <= score_limit:
            break
        yield path_density, path[1:-1]
        # subtract the coverage of this path from the graph and recompute 
        # the edge weights
        subtract_path(G, path[1:-1], path_density)
        # remove path node attributes from graph before
        # calling dynamic programming algorithm again
        clear_path_attributes(G)
    # remove dummy nodes from graph
    G.remove_node(DUMMY_START_NODE)
    G.remove_node(DUMMY_END_NODE)

def assemble_subgraph(G, strand, fraction_major_path, max_paths):
    # find start and end nodes in graph
    tss_id = 0
    start_nodes, end_nodes = find_start_and_end_nodes(G, strand)
    logging.debug("START NODES (TSSs): %s" % start_nodes)
    logging.debug("END NODES: %s" % end_nodes)
    for density, path in find_suboptimal_paths(G, start_nodes, end_nodes,
                                               fraction_major_path=fraction_major_path,
                                               max_paths=max_paths):
        # get tss_id from path
        tss_id = G.node[path[0]][NODE_TSS_ID]
        yield tss_id, 1.0e3 * density, path


def assemble_transcript_graph(G, fraction_major_path, max_paths):
    if fraction_major_path <= 0:
        fraction_major_path = 1e-8        
    # separate graph into connected components
    for Gsub in nx.weakly_connected_component_subgraphs(G):
        # determine proportions of positive and negative strand
        # coverage within nodes in this subgraph
        pos_frac = calculate_strand_fraction(Gsub)
        # transform transcript graph into strand-specific graphs
        # for forward and reverse strands, and calculate node and
        # edge weights
        GG = build_strand_specific_graphs(G, pos_frac)
        

    
    # transform transcript graph into strand-specific graphs
    # for forward and reverse strands, and calculate node and
    # edge weights
    GG = build_strand_specific_graphs(G)
    # assemble strands one at a time
    for strand, Gstrand in enumerate(GG):
        logging.debug("STRAND: %s" % strand)
        gene_id = 0        
        # get connected components - unconnected components are considered
        # different genes
        for Gsub in nx.weakly_connected_component_subgraphs(Gstrand):
            for tss_id, score, path in assemble_subgraph(Gsub, strand, 
                                                         fraction_major_path=fraction_major_path, 
                                                         max_paths=max_paths):
                yield strand, gene_id, tss_id, score, path
            gene_id += 1
