import numpy as np
from numpy import *
import netutils as nu
import itertools as it
import matplotlib.pyplot as plt
import utils.plots as myplots

import netwriter as nw
min_thr = .1
sib_thr = .7
default_name = 'unsup'

def trg_markov_blanket(tkey, tgs = True):
    trgs, tfs = nu.parse_net()
    tg = trgs[tkey]
    
    blnk_tfs = (tg['tfs'],tg['weights'])
    b_tfs = [(blnk_tfs[0][i], blnk_tfs[1][i]) for i in range(len(blnk_tfs[0]))]
    b_tfs = dict(b_tfs)

    if do_tgs:
        blnk_tgs = ([tkey],[1.0])
        for i in range(len(blnk_tfs[0])):
            elt = (blnk_tfs[0][i],blnk_tfs[1][i])
            w0 =elt[1]
            tgs = tfs[elt[0]]['targets']
            ws = array(tfs[elt[0]]['weights'])
            
            blnk_tgs[0].extend(tgs)
            blnk_tgs[1].extend(ws*w0)

        tgn = array(blnk_tgs[0])
        tgw = array(blnk_tgs[1])
        ksrt =argsort(tgn)
        blnk_tgs =[[tgn[i],tgw[i]] for i in ksrt]
    
        t_weights = {}
        for k,g in it.groupby(blnk_tgs,lambda x: x[0]):
            l = list(g)
            ws =array( map(lambda x: x[1],l))
            t_weights[k] = (1 -product(1 - ws))

        b_tgs = t_weights
    else:
        b_tgs = {}
    return b_tgs, b_tfs
        

def coreg_keys(t0 = None, do_plot = False):
    trgs, tfs = nu.parse_net()
    btgs, btfs = trg_markov_blanket(t0)
    if do_plot:
        show_m(btgs,btfs,t0)
    
    min_wt = .3
    tgs_thr =nonzero( greater(btgs.values(),min_wt) )[0]
    keys_thr =[btgs.keys()[i] for i in tgs_thr]
    keys_thr.remove(t0)
    
    na = nu.net_affinity()
    ktf = nu.net_tf_keyidxs()
    ktg = nu.net_trg_keyidxs()

    shared = []
    threshold_sims = True
    for k in keys_thr:
    #    bg,bf = trg_markov_blanket(k, do_tgs = False)#
        row1 = na[ktg[k]] 
        row0 = na[ktg[t0]]
        if threshold_sims:
            row0=array(greater(row0,min_wt),float)
            row1=array(greater(row1,min_wt),float)
        for r in [row0,row1]: 
            l = sqrt(sum(power(r,2)))
            if l != 0: r /= l
        shared.append(sum(row0*row1))

    shared = array(shared)
    min_sharing = .4
    coreg_keys = [ keys_thr[i] for i in nonzero(greater(shared,min_sharing))[0]]
    
    if do_plot:
        plot_shared(shared)

def all_siblings(name = default_name, reset = 0):
    
    donp = True
    hardcopy = True
    if not reset:
        out, sxs = nw.rn2(default_name, 
                          np = donp, hardcopy = hardcopy)
    if reset or not sxs:
        
        trgs, tfs = nu.parse_net(reset=mod(reset,2))
        na = nu.net_affinity(reset = mod(reset,2))
    

    
        na_thr = greater(na, min_thr)
        nrms =  sqrt(sum(power(na_thr,2),1))[:,newaxis]
        nrms[equal(nrms,0)] = 1
        nnn = array(na_thr,float)/nrms
        gg = dot(nnn,nnn.T)

        sibs = array(greater(gg, sib_thr), bool)
    
        nw.wn2(default_name, sibs,
               np = donp,
               hardcopy = hardcopy)
        out = sibs
    return out

def all_siblings_punish_hubs(name = default_name, reset = 0):
    
    donp = True
    hardcopy = True
    if not reset:
        out, sxs = nw.rn2(default_name, 
                          np = donp, hardcopy = hardcopy)
    if reset or not sxs:
        
        trgs, tfs = nu.parse_net(reset=mod(reset,2))
        na = nu.net_affinity(reset = mod(reset,2))

    
        na_thr = greater(na, min_thr)
        
        tfnrms = sqrt(sum(power(na_thr,2),0))[newaxis,:]
        tfnrms[equal(tfnrms,0)] = 1
        na_nohubs = array(na_thr,float) / tfnrms

        nrms =  sqrt(sum(power(na_nohubs,2),1))[:,newaxis]
        nrms[equal(nrms,0)] = 1
        nnn = na_nohubs/nrms
        gg = dot(nnn,nnn.T)
        
        sibs = array(greater(gg, sib_thr), bool)
    
        nw.wn2(default_name, sibs,
               np = donp,
               hardcopy = hardcopy)
        out = sibs
    return out
        
def sib_lists(name = default_name, reset = 0, punish = True):

    if not reset:
        out, sxs = nw.rn2(name = name, hardcopy = True, np = True)
        if not sxs: raise Exception()
    else:
        if punish:
            sibs = all_siblings_punish_hubs(name, reset = mod(reset,2))
        else:
            sibs = all_siblings(name, reset = mod(reset,2))

        counts = zeros(len(sibs))
        for i in sibs:
            counts[i] = len(nonzero(i)[0])
    
        f = plt.figure(2)
        f.clear()
        ax = f.add_subplot(111)
        ax.plot(sort(counts))
        nw.wn2(name , sibs, hardcopy  = True, np = True)
        out = sibs
        
    return out 


def shared_tfs():
    

    raise Exception()
    
    #ck = []
    #ct = 0 
    #for t in trgs:
    #    crk = coreg_keys(t)
    #    ck.append(crk)
    #    ct += 1
    #    if not mod(ct, 10): print ct
    #raise Exception

def plot_shared(shared):
    f = plt.figure(2)
    f.clear()
    ax1 = f.add_subplot(111)
    ax1.plot(sort(shared))

    print len(tgs_thr)
    

    
def show_m(btgs,btfs,name):
    f = plt.figure(1)
    f.clear()
    ax1 = f.add_subplot(211)
    ax1.plot(sorted(btfs.values()))
    ax2 = f.add_subplot(212)
    ax2.plot(sorted(btgs.values()))
    myplots.maketitle(ax1, 'TFS in the Markov blanket for '+name)
    myplots.maketitle(ax2, 'TGS in the Markov blanket for '+name)


def predict(reset = 0):
    trgs, tfs = nu.parse_net(reset = mod(reset,2))
    
    tgk = trgs.keys()
    tfk= tfs.keys()
    ntg, ntf = len(tgk),len(tfk)
    indeg, outdeg = zeros(ntg),zeros(ntf)
    for i in range(ntg): indeg[i] = len(trgs[tgk[i]]['tfs'])
    for i in range(ntf): outdeg[i] = len(tfs[tfk[i]]['targets'])

    
    f = plt.figure(1)
    f.clear()
    ax = f.add_subplot(211)
    ax.plot(indeg[argsort(indeg)])
    ax2 = f.add_subplot(212)
    ax2.plot(outdeg[argsort(outdeg)])

    raise Exception()
