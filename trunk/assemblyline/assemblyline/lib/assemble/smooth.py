'''
Created on Aug 20, 2012

@author: mkiyer
'''
import networkx as nx

from base import NODE_SCORE, SMOOTH_FWD, SMOOTH_REV, SMOOTH_TMP

def smooth_iteration(G, score_attr, smooth_attr):
    nodes = nx.topological_sort(G)
    for u in nodes:
        ud = G.node[u]
        smooth_score = ud[smooth_attr]
        succ = G.successors(u)
        if len(succ) == 0:
            continue
        total_nbr_score = sum(G.node[v][score_attr] for v in succ)        
        if total_nbr_score == 0:
            # if all successors have zero score apply smoothing evenly
            avg_score = smooth_score / len(succ)
            for v in succ:
                vd = G.node[v]
                vd[SMOOTH_TMP] += avg_score
                vd[smooth_attr] += avg_score
        else:
            # apply smoothing proportionately
            for v in succ:
                vd = G.node[v]
                frac = vd[score_attr]/float(total_nbr_score)
                adj_score = frac * smooth_score
                vd[SMOOTH_TMP] += adj_score
                vd[smooth_attr] += adj_score

def smooth_graph(G, score_attr=NODE_SCORE):
    # smooth in forward direction
    smooth_iteration(G, score_attr, SMOOTH_FWD)
    # smooth in reverse direction
    G.reverse(copy=False)
    smooth_iteration(G, score_attr, SMOOTH_REV)
    G.reverse(copy=False)
    # apply densities to nodes
    for n,d in G.nodes_iter(data=True):
        d[score_attr] += d[SMOOTH_TMP]
