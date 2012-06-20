'''
Created on Dec 6, 2010

@author: mkiyer

Algorithm is outlined as following: 
(1) each edge is a 1-clique by itself; 
(2) for k = 2, ..., N, try to expand each (k-1)-clique to k-clique: 
  (2.1) consider a (k-1)-clique the current k-clique KC; 
  (2.2) repeat the following: if for all nodes j in KC, there exists an edge(v,j) add node v to KC
  (2.3) eliminate duplicates; 
(3) the whole graph is N-clique. 
'''

def kcliques(G, k):
    # get adjacency list representation of graph along with 
    # list of (degree, node) tuples
    degree_list = []
    nnbrs = {}
    for n,nbrs in G.adjacency_iter():
        d = len(nbrs)
        degree_list.append((d, n))
        nnbrs[n] = set(nbrs)
    # sort list so nodes with highest degree 
    # appear first
    #deg_sorted_nodes = sorted(degree_list, key=operator.itemgetter(0))
    # each edge is a 1-clique by itself    
    cliques = set(frozenset([n]) for n in nnbrs)
    candidate_cliques = set(frozenset([n]) for n in nnbrs if len(nnbrs[n]) > 0)    
    currentk = 1
    while (len(candidate_cliques) > 0) and (currentk < k):
        rmcliques = set()
        newcliques = set()        
        newcand = set()
        # extend each (k-1)-clique to a k-clique
        for clique in candidate_cliques:
            # intersect neighbors of each node in clique            
            it = iter(clique)
            nbrs = nnbrs[it.next()]
            min_degree = len(nbrs)             
            # subtract nodes that are already in the clique
            cnbrs = set(nbrs).difference(clique)            
            if len(cnbrs) > 0:
                for n in it:
                    nbrs = nnbrs[n]
                    d = len(nbrs)
                    if d < min_degree:
                        min_degree = d
                    cnbrs.intersection_update(nbrs)
                    if len(cnbrs) == 0:
                        break
            # any remaining nodes are part of a larger clique
            if len(cnbrs) > 0:
                # mark existing clique for removal
                rmcliques.add(clique)
                # create larger cliques
                for n in cnbrs:
                    newcliques.add(clique.union([n]))
                # see if this clique is eligible for further joining
                if min_degree > currentk:
                    newcand.update(newcliques)
                #print 'currentk', currentk, 'clique', clique, 'min degree', min_degree
        #print 'candidates', newcand
        cliques.difference_update(rmcliques)
        cliques.update(newcliques)
        # looking for larger cliques
        candidate_cliques = newcand
        currentk += 1      
    return cliques
