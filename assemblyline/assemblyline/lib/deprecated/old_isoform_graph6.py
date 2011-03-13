'''
Created on Feb 20, 2011

@author: mkiyer
'''
def predecessors_iter(G, n, strand):
    for p in G.predecessors(n):
        if p != DUMMY_START_NODE:
            # ignore edges that are not compatible with the
            # strand being considered
            if (strand == POS_STRAND) and (p.end > n.start):
                continue
            elif (strand == NEG_STRAND) and (p.start < n.end):
                continue
        # found a valid predecessor
        yield p

def successors_iter(G, n, strand):
    for s in G.successors(n):
        if s != DUMMY_END_NODE:
            # ignore edges that are not compatible with the
            # strand being considered
            if (strand == POS_STRAND) and (s.start < n.end):
                continue
            elif (strand == NEG_STRAND) and (s.end > n.start):
                continue
        # found a valid successor
        yield s

def get_in_degree(G, n, strand):
    return len(list(predecessors_iter(G, n, strand)))

def get_out_degree(G, n, strand):
    return len(list(successors_iter(G, n, strand)))

def get_strands(G, n):
    return set(x.strand for x in G.node[n]['data'])

def find_start_and_end_nodes(G, strand):
    # find unique starting positions and their 
    # corresponding nodes in the graph
    tss_node_dict = collections.defaultdict(lambda: [])
    end_nodes = []
    for n in G.nodes_iter():
        # only consider nodes on the strand being assembled
        if not all(cmp_strand(strand, x) for x in get_strands(G, n)):
            continue
        # look for start nodes
        in_degree = get_in_degree(G, n, strand)
        if in_degree == 0:
            tss_pos = n.end if strand == NEG_STRAND else n.start
            tss_node_dict[tss_pos].append(n)
            print 'START', 'STRAND', strand, 'NODE', n 
        # look for end nodes
        out_degree = get_out_degree(G, n, strand)
        if out_degree == 0:
            end_nodes.append(n)
            print 'END', 'STRAND', strand, 'NODE', n 
    return tss_node_dict.values(), end_nodes

def calculate_node_weight(exon_data_list):
    # find coverage of node by strand
    strand_scores = [0, 0, 0]
    for edata in exon_data_list:
        strand_scores[edata.strand] += edata.score
    total_score = strand_scores[POS_STRAND] + strand_scores[NEG_STRAND]
    if total_score == 0:
        return strand_scores[NO_STRAND], strand_scores[NO_STRAND]
    # proportionally assign unstranded coverage based on amount of
    # plus and minus strand coverage
    pos_frac = strand_scores[POS_STRAND] / float(total_score)
    pos_weight = strand_scores[POS_STRAND] + pos_frac * strand_scores[NO_STRAND]
    neg_weight = strand_scores[NEG_STRAND] + (1.0 - pos_frac) * strand_scores[NO_STRAND]
    return pos_weight, neg_weight

def calculate_weights(G, strand):
    '''
    calculate edge weights as the fraction of reads-per-million flowing
    out of a given parent node 
    
    edge weights will be stored in the 'weight' attribute of each edge
    '''
    for n,ndict in G.nodes_iter(data=True):
        # find node coverage
        strand_weights = calculate_node_weight(ndict['data'])
        ndict['weight'] = strand_weights[strand]
        # find total coverage flowing out of this node
        total_score = 0
        for succ in successors_iter(G, n, strand):
            intron_data_list = G.edge[n][succ]['data']
            score = sum(x.score for x in intron_data_list)            
            G.edge[n][succ]['score'] = score 
            total_score += score
        # set weight to percentage of coverage flowing out
        for succ in successors_iter(G, n, strand):
            edict = G.edge[n][succ]
            edict['weight'] = edict['score'] / float(total_score)

def add_dummy_start_end_nodes(G, start_nodes, end_nodes):
    # add 'dummy' tss nodes if necessary
    total_weight = sum(G.node[n]['weight'] for n in start_nodes)    
    G.add_node(DUMMY_START_NODE)
    for start_node in start_nodes:        
        weight = G.node[start_node]['weight']
        logging.debug('adding dummy %s -> %s weight=%f total_weight=%f' % (DUMMY_START_NODE, start_node, weight, total_weight))
        G.add_edge(DUMMY_START_NODE, start_node, weight=(weight/total_weight))
    # add a single 'dummy' end node
    G.add_node(DUMMY_END_NODE)
    for end_node in end_nodes:
        logging.debug('adding dummy %s -> %s weight=%f total_weight=%f' % (end_node, DUMMY_END_NODE, total_weight, total_weight))
        G.add_edge(end_node, DUMMY_END_NODE, weight=1.0)

# constants for use as graph attributes
PLEN = 'length'
PCOV = 'weight'
PFRAC = 'frac'
PSCORE = 'score'
PSRC = 'source'
EDGE_WEIGHT = 'weight'

def visit_node_from_parent(G, parent, child):    
    # compute length, coverage, and rpkm of the current
    # path when extended to include the child
    pattrs = G.node[parent]
    # new path length is parent length + child length    
    path_length = pattrs.get(PLEN, 0) + (child.end - child.start)    
    # new path coverage is parent coverage + child coverage
    path_cov = pattrs.get(PCOV, 0) + G.node[child].get(PCOV, 0)
    # keep track of the minimum fraction of coverage that
    # flows between and two nodes along the entire path     
    parent_frac = pattrs.get(PFRAC, 1.0)
    child_frac = G.edge[parent][child][EDGE_WEIGHT]
    path_frac = parent_frac if parent_frac < child_frac else child_frac
    # score is total coverage multiplied by the path fraction
    # divided by the length of the path
    path_score = (path_frac * path_cov) / path_length
    return path_length, path_cov, path_frac, path_score

def dyn_prog_search(G, source, strand):
    """Find the highest scoring path by dynamic programming"""
    # Taken from NetworkX source code http://networkx.lanl.gov
    stack = [source]
    while stack:
        parent = stack[0]
        for child in successors_iter(G, parent, strand):            
            logging.debug('Visiting %s -> %s' % (parent, child))            
            path_length, path_cov, path_frac, path_score = visit_node_from_parent(G, parent, child)
            logging.debug('\tlen=%d frac=%f score=%f child_score=%f' %
                          (path_length, path_frac, path_score, G.node[child].get(PSCORE,-1)))
            # only update this node if it has not been visited
            # or if the new score is higher than the old score
            cattrs = G.node[child]
            if path_score > cattrs.get(PSCORE, -1):
                # keep pointer to parent node that produced this high scoring path
                cattrs[PLEN] = path_length
                cattrs[PCOV] = path_cov
                cattrs[PFRAC] = path_frac
                cattrs[PSCORE] = path_score
                cattrs[PSRC] = parent
                logging.debug('\tupdated child len=%d cov=%f frac=%f score=%f' % 
                              (path_length, path_cov, path_frac, path_score))
                # continue iterating through the child path
                stack.append(child)
        stack.pop(0)

def find_best_path(G, source, sink, strand):
    '''
    use dynamic programming to find the highest scoring path through 
    the graph starting from 'source'
    '''
    # run the dynamic programming algorithm to
    # score paths
    dyn_prog_search(G, source, strand)
    # traceback
    path = [sink]
    pattrs = G.node[sink]
    score = pattrs[PSCORE]
    frac = pattrs[PFRAC]    
    length = pattrs[PLEN]
    cov = pattrs[PCOV]
    while path[-1] != source:
        #logging.debug('path=%s node=%s attrs=%s' %
        #              (path, path[-1], G.node[path[-1]]))
        path.append(G.node[path[-1]][PSRC])
    path.reverse()
    #logging.debug("FINAL path=%s" % (path))
    # clear path attributes
    for n in G.nodes_iter():
        nattrs = G.node[n]
        if PLEN in nattrs:
            del nattrs[PLEN]
            del nattrs[PCOV]
            del nattrs[PFRAC]
            del nattrs[PSCORE]
            del nattrs[PSRC]
    return path, score, frac, cov, length


def find_suboptimal_paths(G, start_node, end_node, strand, 
                          fraction_major_path,
                          max_paths,
                          max_iters=10000):
    # copy the EDGEWT attribute rather than the entire graph
    edge_weight_dict = {}
    for u,v,d in G.edges_iter(data=True):
        if EDGE_WEIGHT in d:
            edge_weight_dict[(u,v)] = d[EDGE_WEIGHT]
    # find overall highest scoring path through graph
    path, path_score, path_frac, path_cov, path_length = \
        find_best_path(G, start_node, end_node, strand)
    covscore = 1.0e3 * path_cov * path_frac / path_length
    score_limit = covscore * fraction_major_path
    logging.debug("Maximum path covscore=%f limit=%f weight=%f psiscore=%f path=%s" % 
                  (covscore, score_limit, path_frac, path_score, path))
    # enumerate paths until the score falls below the 
    # specified percentage of the best path
    while covscore > score_limit:
        yield covscore, path
        # subtract path weight from edge weights
        for parent,child in zip(path[:-1], path[1:]):
            pcedge = G.edge[parent][child]
            if EDGE_WEIGHT in pcedge:
                pcedge[EDGE_WEIGHT] -= path_frac
        # get next path        
        path, path_score, path_frac, path_cov, path_length = \
            find_best_path(G, start_node, end_node, strand)
        covscore = 1.0e3 * path_cov * path_frac / path_length
        #logging.debug("Current path covscore=%f limit=%f frac=%f score=%f path=%s" % 
        #              (covscore, score_limit, path_frac, path_score, path))
        logging.debug("Current path covscore=%f limit=%f frac=%f score=%f" % 
                      (covscore, score_limit, path_frac, path_score))
    # restore original edge weights
    for u,v,d in G.edges_iter(data=True):
        if EDGE_WEIGHT in d:
            d[EDGE_WEIGHT] = edge_weight_dict[(u,v)]

def assemble_strand(G, strand,
                    fraction_major_isoform,
                    max_paths,
                    max_iters):
    # calculate weights for edges
    calculate_weights(G, strand)
    # find start and end nodes in graph
    tss_id = 0
    start_node_lists, end_nodes = find_start_and_end_nodes(G, strand)
    for start_nodes in start_node_lists:
        logging.debug("STRAND: %s" % strand)
        logging.debug("START NODES (TSSs): %s" % start_nodes)
        logging.debug("END NODES: %s" % end_nodes)
        # at an artificial node at the start and end that all paths with
        # the same tss are searched together
        add_dummy_start_end_nodes(G, start_nodes, end_nodes)
        for score, path in find_suboptimal_paths(G, DUMMY_START_NODE, 
                                                 DUMMY_END_NODE,
                                                 strand,
                                                 fraction_major_path=fraction_major_isoform,
                                                 max_paths=max_paths,
                                                 max_iters=max_iters):
            # remove dummy nodes (first and last) when returning path
            yield tss_id, score, path
        #path[1:-1]
        # remove dummy nodes from graph
        G.remove_node(DUMMY_START_NODE)
        G.remove_node(DUMMY_END_NODE)
        tss_id += 1