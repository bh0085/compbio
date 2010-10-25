import netsvd
import matplotlib.pyplot as plt
from numpy import *
import numpy as np
import utils.colors as mycolors
import netutils as nu

def draw_genegene(gg):
    maxview = 200
    f = plt.figure(0)
    ax = f.add_axes([0,0,1,1])
    ggd = matrix(gg[0:maxview,0:maxview])
    #for i in range(len(ggd)):
    #    ggd[i,i] = 0

    ax.imshow(ggd)

def draw_hclustered(clustered):
    
    c = clustered['HCluster']
    clusters = c.cut(0)
    f = plt.figure(0)
    f.clear()
    ax = f.add_axes([0,0,1,1],aspect = 'auto')
    ax.set_aspect('auto')
    c0 = c.cut(0)
    n = max(array(c0)) +1
    nlvl = 100
    h = zeros((nlvl,n))

    cuts = linspace(0,.9,nlvl)
    appearances = zeros(n)-1
    hinds = zeros((nlvl,n),int)
    for i in range(nlvl):
        cut = cuts[i]
        clusters = c.cut(cut)
        for j in range(n):
            cval = clusters[j]
            if appearances[cval] == -1: appearances[cval] = j
            h[i,j] = appearances[cval]
        hinds[i,:] = argsort(h[i,:])
        h[i,:]=h[i,hinds[i,:]]

    for i in range(shape(h)[0]):
        h[i,:] /= max(h[i,:])

    ax.imshow(mycolors.imgcolor(h, BW = True), aspect = 'auto', interpolation = 'nearest')
    saff,ks = nu.net_square_affinity()
    raff,kt,ktf = nu.net_affinity()
    
    tfidxs = []
    for tf in ktf.keys():
        tfidxs.append(ks.index(tf))
    tfidxs= array(tfidxs)

    is_clustered = tfidxs[nonzero(less(tfidxs,n))[0]]
    ntf = len(is_clustered)
    tf_alpha = zeros((nlvl,n))
    for i in range(ntf):
        tf_alpha +=equal(hinds,is_clustered[i])
        
    tf_rgba = mycolors.imgcolor(tf_alpha,alpha = True,color = [1,0,0])
    
    ax.imshow(tf_rgba, aspect = 'auto', interpolation = 'nearest')

