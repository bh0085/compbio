import pickle
from numpy import *
import matplotlib.pyplot as plt
import netutils as nu
import numpy as np
import numpy.linalg as lin
def compute(namedefault_name, number = 0, reset = False):
    #possible names: bRN, kRN, fRN, mRN
    N = nu.net_affinity(name, reset = reset)
    U, S, Vh = lin.svd(N)
    V = Vh.T
    uvecs = U[:,:len(S)]
    vvecs = V
    dosave = ((kmap0,kmap1),(uvecs, S, vvecs))
    
    savename = name + '_svd.pickle'
    return savename

def get_SVD(name, reset = False):
    if reset: 
        savename = compute(name = name, reset = reset) 
    else:
        savename = name + '_svd.pickle'
    return nu.readnet(savename)

def view_net(name = default_name):
    keys,svd = get_SVD(name,reset = False)
    U, S, V = svd
    view_svd(U,S,V, fig = i)

def view_svd(U,S,V, fig = 0 ):    
    m = shape(U)[0]
    n = shape(V)[0]

    nt = shape(U)[0]
    yvals = pow(reshape(mod(arange(m*n) , n),(m,n)),2)
    tots = sum( abs(U) *yvals,1)
    srt = argsort(tots)
    Usrt = U[srt,:]
        
    fig = plt.figure(fig)
    fig.clear()
    ax = fig.add_axes([0,0,1,1],frameon = False)
    ax.imshow(abs(U[srt[0:m/2],:]),aspect = 'auto',interpolation = 'nearest')

