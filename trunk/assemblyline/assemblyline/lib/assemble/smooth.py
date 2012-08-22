'''
Created on Aug 20, 2012

@author: mkiyer
'''
import networkx as nx

from base import NODE_DENSITY, SMOOTH_FWD, SMOOTH_REV, SMOOTH_TMP

def smooth_iteration(G, density_attr, smooth_attr):
    nodes = nx.topological_sort(G)
    for u in nodes:
        ud = G.node[u]
        smooth_density = ud[smooth_attr]
        succ = G.successors(u)
        if len(succ) == 0:
            continue
        total_nbr_density = sum(G.node[v][density_attr] for v in succ)        
        if total_nbr_density == 0:
            # if all successors have zero density apply smoothing evenly
            avg_density = smooth_density / len(succ)
            for v in succ:
                vd = G.node[v]
                vd[SMOOTH_TMP] += avg_density
                vd[smooth_attr] += avg_density
        else:
            # apply smoothing proportionately
            for v in succ:
                vd = G.node[v]
                frac = vd[density_attr]/float(total_nbr_density)
                adj_density = frac * smooth_density
                vd[SMOOTH_TMP] += adj_density
                vd[smooth_attr] += adj_density

def smooth_graph(G, density_attr=NODE_DENSITY):
    # smooth in forward direction
    smooth_iteration(G, density_attr, SMOOTH_FWD)
    # smooth in reverse direction
    G.reverse(copy=False)
    smooth_iteration(G, density_attr, SMOOTH_REV)
    G.reverse(copy=False)
    # apply densities to nodes
    for n,d in G.nodes_iter(data=True):
        d[density_attr] += d[SMOOTH_TMP]
