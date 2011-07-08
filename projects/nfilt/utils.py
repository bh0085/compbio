import networkx as nx
from numpy import *
import numpy as np
import scipy.linalg as slg
import cb.utils.memo as mem

import cb.utils.plots as myplots

def nx_from_network(targets, nodes_allowed = None,
                    name = 'default', **kwargs):
    def set_nx(**kwargs):
        targets = kwargs.get('targets')
        nodes_allowed = kwargs.get('nodes_allowed')
 
        if nodes_allowed != None:
            edges = [(tf, tg, float(wt)) 
                     for tg, elt in targets.iteritems() if tg in nodes_allowed
                     for tf, wt  in zip(elt['tfs'], elt['weights']) if tf in nodes_allowed]
        else:
            edges = [(tf, tg, float(wt)) 
                     for tg, elt in targets.iteritems()
                     for tf, wt  in zip(elt['tfs'], elt['weights']) ]
    
        graph = nx.DiGraph()
        graph.add_weighted_edges_from(edges)
        return graph

    return mem.getOrSet(set_nx,
                        **mem.rc(kwargs,
                                 register = name,
                                 targets = targets, 
                                 nodes_allowed = nodes_allowed, 
                                 on_fail = 'compute'))

def dotprod(adj1, adj2):
    return sum(adj1 * adj2)


def cosine_adj(adj1,adj2):

    cosine = sum(adj1*adj2)
    if cosine > 1.01:
        raise Exception()
    return sum(adj1*adj2)

def thr_graph(g1, thr = .9):
    g2 = type(g1)()
    g2.add_nodes_from(g1.nodes())
    g2.add_weighted_edges_from([ (e[0],e[1],1.)
                        for e in nx.to_edgelist(g1) 
                                 if e[2]['weight'] >= thr])
    
    return g2

def top_edges(g1, max_edges = 100):
    g2 = type(g1)()
    g2.add_nodes_from(g1.nodes())
    
    ewts = [x[2]['weight'] for x in nx.to_edgelist(g1)]
    keepwts = set(argsort(ewts)[::-1][:max_edges])
    g2.add_weighted_edges_from([ (e[0],e[1],e[2]['weight'])
                                 for i,e in enumerate(nx.to_edgelist(g1))
                                 if i in keepwts])
    
    return g2    

def filter_sparse(g1, n_c = 5, 
                  max_edges = -1,
                  last_component = False):
    
    '''
    Filter a sparse version of the network by PCA.

g1:  The input network graph.
n_c: The number of principal components to compute
max_edges: The maximum number of edges to keep. -1 => keep all.
last_component: Keep only the final principal component
'''

    import scipy.sparse.linalg as las
    import scipy.sparse.lil as ll
    import scipy.sparse as ssp

    adj = ssp.csr_matrix(nx.to_scipy_sparse_matrix(g1))
    nodes = g1.nodes()

    U,s, Vh = svd = las.svd(adj, n_c)
    

    U[less(abs(U), .001)] = 0
    Vh[less(abs(Vh), .001)] = 0
    
    if last_component:
        s_last = s; s_last[1:] *= 0
        filtered = ll.lil_matrix(U)*ll.lil_matrix(diag(s_last)) *ll.lil_matrix(Vh)
    else:
        filtered = ll.lil_matrix(U)*ll.lil_matrix(diag(s)) *ll.lil_matrix(Vh)

    if max_edges != -1:
        filtered.data[argsort(abs(filtered.data))[:-1 * max_edges]] = 0
        filtered.eliminate_zeros()


    g = nx.DiGraph() 
    g.add_nodes_from(nodes)
    g.add_weighted_edges_from([(nodes[nz[0]],nodes[nz[1]],nz[2])
                               for nz in zip(*ssp.find(filtered))])
    
    
    

    
    return g


def filter_graph(g1, n_c = 5):
    adj = nx.adj_matrix(g1)
    nodes = g1.nodes()
    
    U,s,Vh = slg.svd(adj)
    s[n_c:] = 0.
    filtered =  dot(dot(U ,diag(s)),Vh)
    filtered[less(filtered, 1e-8)] = 0
    
    g = nx.DiGraph() 
    g.add_nodes_from(nodes)
    for nz in zip(*nonzero(filtered)):
        g.add_weighted_edges_from([(nodes[nz[0]],nodes[nz[1]],filtered[nz])])
        
    return g
