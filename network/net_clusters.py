
import netsvd as ns
import netutils as nu
import netwriter as nw
import mlpy
import compbio.learning.clusterview as cl
from numpy import *

name = nu.default_name
maxgenes = 8316
k = 10

def runall(maxg = maxgenes, runk = 40, rname = nu.default_name):
    maxgenes = maxg
    k = runk
    name = rname
    

    print '''
...netclusters.runall running with:

k='''+str(k)+'''
maxgenes='''+str(maxgenes)+'''

'''
    r = 2
    net_cluster_gg(k,reset = r)
    net_cluster_ggn(k,reset = r)
    net_cluster_genes_by_tf(k,reset = r)
   
    dosvd = False
    if dosvd:

        net_cluster_ggsvdU(k, reset = r)
    #for now, we are not working witht he clusters of the V arrays
    #I don't understand what they represent.
    #(And while we are clustering V.T, I suspect that we should
    #be clustering V instead. who knows?)
    #net_cluster_ggsvdV(k, reset = r)
        net_cluster_sqsvdU(k, reset = r)
    #net_cluster_sqsvdV(k, reset = r)

def viewall(projs,idxs = [0,1],axis = 'tf',viewmany = False):
    clusters = []
    means = []
    kmeans = []
    kmeans.append(net_cluster_gg())
    kmeans.append(net_cluster_ggn())
    kmeans.append(net_cluster_genes_by_tf())
    dosvd = False
    if dosvd:
        kmeans.append(net_cluster_ggsvdU())
    #kmeans.append(net_cluster_ggsvdV())
        kmeans.append(net_cluster_sqsvdU())
    #kmeans.append(net_cluster_sqsvdV())
    
    clusters = map(lambda x: x[0], kmeans)
    means = map(lambda x: x[1], kmeans)

    if viewmany:
        projs = cl.viewmany(means,clusters,fig = 13)
    else:
        cl.one(means,clusters,projs,idxs = idxs, axis = axis)
    
    return projs


def net_cluster_gg(k = k, reset = 0):
    hardcopy = True
    try:
        if reset: raise Exception('compute')
        out,sxs =  nw.rn2(name, hardcopy = hardcopy)
        if not sxs: raise Exception()

    except Exception as e:
        if e.args[0] != 'compute': raise Exception()
        nw.claim_reset()    
        gg = nu.net_genegene(reset = mod(reset,2))
        kmeans = mlpy.Kmeans(k)
        clustered = kmeans.compute(gg[0:maxgenes,:])
        means = kmeans.means
        out = (clustered,means)
        nw.wn2(name, (clustered, means) ,hardcopy = hardcopy)

    return out
def net_cluster_ggn(k = k, reset = 0):
    hardcopy = True
    try:
        if reset: raise Exception('compute')
        out,sxs =  nw.rn2(name, hardcopy = hardcopy)
        if not sxs: raise Exception()

    except Exception as e:
        if e.args[0] != 'compute': raise Exception()
        nw.claim_reset()   
        gg = nu.net_genegene_norm(reset = mod(reset,2))
        kmeans = mlpy.Kmeans(k)
        clustered = kmeans.compute(gg[0:maxgenes,:])
        means = kmeans.means
        out = (clustered,means)
        nw.wn2(name, (clustered, means) ,hardcopy = hardcopy) 
    return out

def net_cluster_genes_by_tf(k= k, reset = 0):
    hardcopy = True
    try:
        if reset: raise Exception('compute')
        out,sxs =  nw.rn2(name, hardcopy = hardcopy)
        if not sxs: raise Exception()

    except Exception as e:
        if e.args[0] != 'compute': raise Exception()
        nw.claim_reset()    
        gg = nu.net_affinity(reset = mod(reset,2))
        kmeans = mlpy.Kmeans(k)
        clustered = kmeans.compute(gg[0:maxgenes,:])
        means = kmeans.means
        out = (clustered,means)
        nw.wn2(name, (clustered, means) ,hardcopy = hardcopy) 
    return out

def net_cluster_ggsvdU(k = k,reset = 0):
    hardcopy = True
    try:
        if reset: raise Exception('compute')
        out,sxs =  nw.rn2(name, hardcopy = hardcopy)
        if not sxs: raise Exception()


    except Exception as e:
        if e.args[0] != 'compute': raise Exception()
        nw.claim_reset() 

        gg = ns.net_ggsvd_U()
        kmeans = mlpy.Kmeans(k)
        clustered = kmeans.compute(gg[0:maxgenes,:])
        means = kmeans.means
        nw.wn2(name, (clustered, means) ,hardcopy = hardcopy) 
    return out

def net_cluster_ggsvdV(k = k,reset = 0):
    hardcopy = True
    try:
        if reset: raise Exception('compute')
        out,sxs =  nw.rn2(name, hardcopy = hardcopy)
        if not sxs: raise Exception()

    except Exception as e:
        if e.args[0] != 'compute': raise Exception()
        nw.claim_reset()
        gg = ns.net_ggsvd_V().T
        kmeans = mlpy.Kmeans(k)
        clustered = kmeans.compute(gg[0:maxgenes,:])
        means = kmeans.means
        nw.wn2(name, (clustered, means) ,hardcopy = hardcopy) 
    return out

def net_cluster_sqsvdU(k = k,reset = 0):
    hardcopy = True
    try:
        if reset: raise Exception('compute')
        out,sxs =  nw.rn2(name, hardcopy = hardcopy)
        if not sxs: raise Exception()



    except Exception as e:
        if e.args[0] != 'compute': raise Exception()
        nw.claim_reset() 

        gg = ns.net_square_svd_U()
        kmeans = mlpy.Kmeans(k)
        clustered = kmeans.compute(gg[0:maxgenes,:])
        means = kmeans.means
        nw.wn2(name, (clustered, means) ,hardcopy = hardcopy) 
    return out
   
def net_cluster_sqsvdV(k = k,reset = 0):
    hardcopy = True
    try:

        if reset: raise Exception('compute')
        out,sxs =  nw.rn2(name, hardcopy = hardcopy)
        if not sxs: raise Exception()


    except Exception as e:
        if e.args[0] != 'compute': raise Exception()
        nw.claim_reset()
        gg = ns.net_square_svd_V().T
        kmeans = mlpy.Kmeans(k)
        clustered = kmeans.compute(gg[0:maxgenes,:])
        means = kmeans.means
        nw.wn2(name, (clustered, means) ,hardcopy = hardcopy) 
    return out
