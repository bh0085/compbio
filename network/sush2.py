import os
import numpy as np
import re
import itertools as it
import scipy as sp
import scipy.linalg as lin
import matplotlib.pyplot as plt

def parse_net(name = 'mRN',number = 0):
    prefix = os.path.join('predmodel/regressionwts', name)
    nwdata = open(os.path.join(prefix,
                               'nw_'+str(number) + '.sif')).read()
    #A few functions defined here to be used later
    trgfun = lambda x: x[1]
    wtfun = lambda x:float( x[2] )
    tffun = lambda x: x[0]
    sigmafun = lambda x: 1 / (1 + np.exp(-x /1))

    r = re.compile('^[ ]*(?P<target>[^\s]+)\t(?P<tf>[^\s]+)\t(?P<weight>[^\s]+)'
                   ,re.M)
    matches = list(re.finditer(r,nwdata))    
    #Unsorted lists of tfs and targets
    targets =map(lambda x:x.group('target'),matches)
    tfs =    map(lambda x:x.group('tf'),matches)
    weights =map(lambda x:x.group('weight'),matches)
    
    #Concat the data for easier sorting
    cat = []
    for i in np.argsort(tfs):
        cat.append([tfs[i],targets[i],weights[i]])

    #Extract a dictionary with information for each target.
    trg_d = {}
    count = 0.0
    for k, g in it.groupby(sorted(cat,key = trgfun),key = trgfun):
        l = list(g)
        count += 1.0
        trg_d[k] = {'color': np.array([count, 0, 0]),
                    'tfs' : map(tffun,l),
                    'weights': map(wtfun,l)
                    }


    #Extract a dictionary with information for each TF
    tf_d = {}
    for k, g in it.groupby(cat,key = lambda x: x[0]):
        l = list(g)
        tf_targets = map(lambda x: x[1],l)
        
        tf_d[k] = {'targets':tf_targets}

    return (trg_d, tf_d)


def svd_net(name = 'mRN', number = 0):
    #possible names: bRN, kRN, fRN, mRN
    trgs, tfs = parse_net(name, number)
    print 'N Targets: ' + str(len(trgs.keys()))
    print 'N TFs: ' + str(len(tfs.keys()))
    
    k0 = trgs.keys()
    k1 = tfs.keys()

    kmap0 = {}
    kmap1 = {}
    for i in range(len(k0)):
        kmap0[k0[i]] = i
    for i in range(len(k1)):
        kmap1[k1[i]] = i

    #Now, for a first shot, lets construct an n*m matrix
    #with possible interactions for each gene appearing in 
    #the network

    m = len(kmap0.keys())
    n = len(kmap1.keys())
    N = np.zeros([m,n])
    for k_trg,v_trg in trgs.items():
        i = kmap0[k_trg]
        for k_tf in v_trg['tfs']:            
            j = kmap1[k_tf]
            N[i][j] = 1.0

    import pickle
    svdname = name + '.svd'
    comp_svd = True
    if comp_svd:
        U, S, Vh = lin.svd(N)
        V = Vh.T
        svd = (U,S,V)
        pickle.dump(svd,open(svdname,'w'))
    else:
        U,S,V = pickle.load(open(svdname))

    #from mlabwrap import mlab
    #U2, S2, V2 = mlab.svd(N,nout = 3)
    
    f = plt.figure(0)
    plt.clf()
    ax = plt.axes(frameon = False)
    ax.set_xlim([0,1])
    ax.set_ylim([0,1])

    xax = np.arange(float(len(S)))/len(S)
    yax = S / S.max()
    ax.plot(xax, yax)


    print N
