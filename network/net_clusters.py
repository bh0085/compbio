
import netsvd as ns
import netutils as nu
import netwriter as nw
import mlpy
import compbio.learning.clusterview as cl

name = nu.default_name
maxgenes = 8321

def runall(maxg = 8321, runk = 40, rname = nu.default_name):
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
    net_cluster_ggsvdU(k, reset = r)
    #for now, we are not working witht he clusters of the V arrays
    #I don't understand what they represent.
    #(And while we are clustering V.T, I suspect that we should
    #be clustering V instead. who knows?)
    #net_cluster_ggsvdV(k, reset = r)
    net_cluster_sqsvdU(k, reset = r)
    #net_cluster_sqsvdV(k, reset = r)

def viewall(projs,idxs = [0,1],axis = 'tf'):
    clusters = []
    means = []
    kmeans = []
    kmeans.append(net_cluster_gg())
    kmeans.append(net_cluster_ggn())
    kmeans.append(net_cluster_genes_by_tf())
    kmeans.append(net_cluster_ggsvdU())
    #kmeans.append(net_cluster_ggsvdV())
    kmeans.append(net_cluster_sqsvdU())
    #kmeans.append(net_cluster_sqsvdV())
    
    clusters = map(lambda x: x[0], kmeans)
    means = map(lambda x: x[1], kmeans)

    viewmany = False
    if viewmany:
        projs = cl.viewmany(means,clusters,fig = 13)
    else:
        cl.one(means,clusters,projs,idxs = idxs, axis = axis)
    
    return projs


def net_cluster_gg(k = k, reset = 0):
    hardcopy = True
    try:
        if reset: raise Exception('compute')
        return nw.readnet(name, hardcopy = hardcopy)

    except Exception as e:
        if e.args[0] != 'compute': raise Exception()
        nw.claim_reset()    
        gg = nu.net_genegene()
        kmeans = mlpy.Kmeans(k)
        clustered = kmeans.compute(gg[0:maxgenes,:])
        means = kmeans.means
        nw.writenet(name, (clustered, means) ,hardcopy = hardcopy)


def net_cluster_ggn(k = k, reset = 0):
    hardcopy = True
    try:
        if reset: raise Exception('compute')
        return nw.readnet(name, hardcopy = hardcopy)

    except Exception as e:
        if e.args[0] != 'compute': raise Exception()
        nw.claim_reset()   
        gg = nu.net_genegene_norm()
        kmeans = mlpy.Kmeans(k)
        clustered = kmeans.compute(gg[0:maxgenes,:])
        means = kmeans.means
        nw.writenet(name, (clustered, means) ,hardcopy = hardcopy) 

def net_cluster_genes_by_tf(k= k, reset = 0):
    hardcopy = True
    try:
        if reset: raise Exception('compute')
        return nw.readnet(name, hardcopy = hardcopy)

    except Exception as e:
        if e.args[0] != 'compute': raise Exception()
        nw.claim_reset()    
        gg = nu.net_affinity()[0]
        kmeans = mlpy.Kmeans(k)
        clustered = kmeans.compute(gg[0:maxgenes,:])
        means = kmeans.means
        nw.writenet(name, (clustered, means) ,hardcopy = hardcopy) 

def net_cluster_ggsvdU(k = k,reset = 0):
    hardcopy = True
    try:
        if reset: raise Exception('compute')
        return nw.readnet(name, hardcopy = hardcopy)

    except Exception as e:
        if e.args[0] != 'compute': raise Exception()
        nw.claim_reset() 

        gg = ns.net_ggsvd_U()
        kmeans = mlpy.Kmeans(k)
        clustered = kmeans.compute(gg[0:maxgenes,:])
        means = kmeans.means
        nw.writenet(name, (clustered, means) ,hardcopy = hardcopy) 
   
def net_cluster_ggsvdV(k = k,reset = 0):
    hardcopy = True
    try:
        if reset: raise Exception('compute')
        return nw.readnet(name, hardcopy = hardcopy)

    except Exception as e:
        if e.args[0] != 'compute': raise Exception()
        nw.claim_reset()
        gg = ns.net_ggsvd_V().T
        kmeans = mlpy.Kmeans(k)
        clustered = kmeans.compute(gg[0:maxgenes,:])
        means = kmeans.means
        nw.writenet(name, (clustered, means) ,hardcopy = hardcopy) 

def net_cluster_sqsvdU(k = k,reset = 0):
    hardcopy = True
    try:
        if reset: raise Exception('compute')
        return nw.readnet(name, hardcopy = hardcopy)

    except Exception as e:
        if e.args[0] != 'compute': raise Exception()
        nw.claim_reset() 

        gg = ns.net_square_svd_U()
        kmeans = mlpy.Kmeans(k)
        clustered = kmeans.compute(gg[0:maxgenes,:])
        means = kmeans.means
        nw.writenet(name, (clustered, means) ,hardcopy = hardcopy) 
   
def net_cluster_sqsvdV(k = k,reset = 0):
    hardcopy = True
    try:
        if reset: raise Exception('compute')
        return nw.readnet(name, hardcopy = hardcopy)

    except Exception as e:
        if e.args[0] != 'compute': raise Exception()
        nw.claim_reset()
        gg = ns.net_square_svd_V().T
        kmeans = mlpy.Kmeans(k)
        clustered = kmeans.compute(gg[0:maxgenes,:])
        means = kmeans.means
        nw.writenet(name, (clustered, means) ,hardcopy = hardcopy) 
