import cb.p.network.io as nio
import comparator as comp

from numpy import *

import scipy.sparse.linalg as las
import scipy.sparse.lil as ll
import scipy.sparse as ssp
import scipy.linalg as la

import networkx as nx

import cb.utils.plots as myplots

def cmp2(graph):
    '''Run with a graph such as from comparator.get_graphs()['bn']'''
    adj = ssp.csr_matrix(nx.to_scipy_sparse_matrix(graph))
    tmp = array(adj[:1000,:1000].todense())
    
    tf = myplots.tmpF(3,(12,12))
    f = tf.f
    ax = f.add_subplot(111)

    gg = dot(tmp.T, tmp)
    
    eigs = la.eig(gg)
    vecs = eigs[1]
    #raise Exception(imsho)
    ax.imshow(array(eigs[1], float))
    tf.save()

