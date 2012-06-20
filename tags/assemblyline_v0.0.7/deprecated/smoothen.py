'''
Created on Dec 15, 2011

@author: mkiyer
'''
import numpy as np

from assembler_base import NODE_DENSITY, EDGE_OUT_FRAC

def smoothen_graph_leastsq(G):
    # cannot smooth a graph with just one node
    if len(G) == 1:
        return
    # label nodes
    cur_id = 0
    eqns = {}
    cols = {}
    y = np.zeros(len(G), dtype=np.float)
    # assign numeric ids to nodes
    for n,d in G.nodes_iter(data=True):
        d['id'] = cur_id
        y[cur_id] = d[NODE_DENSITY]
        if G.in_degree(n) > 0:
            a = np.zeros(len(G), dtype=np.float)
            a[cur_id] = -1.0
            eqns[cur_id] = a
        cols[cur_id] = n
        cur_id += 1
    # setup equations for each node
    for u,v,d in G.edges_iter(data=True):
        u_id = G.node[u]['id']
        v_id = G.node[v]['id']
        out_frac = d[EDGE_OUT_FRAC]
        #frac = d[EDGE_OUT_FRAC] * d[EDGE_IN_FRAC]
        eqns[v_id][u_id] = out_frac
    # build matrix for least squares error
    # minimization
    row_ids = sorted(eqns)
    A = np.vstack(eqns[id] for id in row_ids)
    A = np.asmatrix(A)
    I = np.identity(len(G), dtype=np.float)
    # perform least squares minimization
    beta = (I - A.T*((A*A.T).I)*A).dot(y)
    # set new node weights
    for col_id,n in cols.iteritems():
        G.node[n][NODE_DENSITY] = beta[0,col_id]
        del G.node[n]['id']
    #for n,d in G.nodes_iter(data=True):
    #    print 'NEW NODE', n, d