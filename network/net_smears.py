import netutils as nu
import netsvd as ns
from numpy import *
import numpy as np
import matplotlib.pyplot as plt
import learning.mycluster as mc
import mlpy

xskip = 8000/4000
yskip = 8000/300
#set an xmax if you are smearing - for eg, gene in degree
xmax = 200
cfig = 0
def smear_sqa(smear_gene = True):
    sqa = nu.net_square_affinity()
    
    if smear_gene:
        img = array(sqa[0].T)
    else:
        img = array(sqa[0])

    s = shape(img)
    for i in range(s[1]):
        img[:,i] = np.sort(img[:,i])[::-1]

    sums = np.sum(img,0)
    img = img[:,np.argsort(sums)[::-1]]
    draw_smear(img[:xmax:xskip,::yskip])    


def cluster( columns = False):

    sqa = nu.net_affinity()[0]
    #tfidxs = nonzero(np.max(sqa ,0))[0]
    #raise Exception()
    
    if columns:
        arr = array(sqa.T)
    else:
        arr = array(sqa)
        
    subarr = arr[:,:]
    
    k = 20
    
    kmeans = mlpy.Kmeans(k)
    clustered = kmeans.compute(subarr)
    means = kmeans.means
    import compbio.learning.clusterview as cl
    cl.specview(means,clustered)

#mc.test_withdata(100,subarr)
    
        
   
def smear_ggn():
    ggn = nu.net_genegene_norm()
    img = array(ggn)
    s = shape(img)
    for i in range(s[1]):
        img[:,i] = np.sort(img[:,i])[::-1]/ np.max(img[:,i])
    sums = np.sum(img ** 2,0)
    img = img[:, np.argsort(sums)[::-1]]
    draw_smear(img[::xskip,::yskip])    


    
def smear_gg():
    ggn = nu.net_genegene()
    img = array(ggn)
    s = shape(img)
    for i in range(s[1]):
        img[:,i] = np.sort(img[:,i])[::-1]/ np.max(img[:,i])
    sums = np.sum(img ** 2,0)
    img = img[:, np.argsort(sums)[::-1]]
    draw_smear(img[::xskip,::yskip])    

def draw_ggn():
    ggn = nu.net_genegene_norm()
    f = plt.figure(cfig)
    f.clear()
    ax = f.add_axes([0,0,1,1])
    sums = np.sum(ggn,0)
    ax.plot(sums)

def draw_smear(smear ):
    f = plt.figure(cfig)
    ax = f.add_axes([0,0,1,1])
    ax.imshow(smear,aspect = 'auto')
    
