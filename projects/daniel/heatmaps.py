import matplotlib.pyplot as plt
from numpy import *
import numpy as np
import compbio.config as cfg

def load(net = 2, num = 1676,
         min_module_size = 10,
         min_go_size = 5,
         max_go_modules = 3,
         prb_threshold = [.01, .5]):
    fopen = open(cfg.dataPath('daniel/heatplots/net{0}_top11}_heatplot_matrix.txt'.format(net,num)),'r')
    l0 = fopen.readline()

    arr = array([[float(elt) for elt in l.split('\t')] for  l in  fopen.xreadlines() if l.strip() != ''])
    ids = arr[:,:1]
    arr = arr[:,1:]

    arr[equal(arr,0)] = np.min(arr[not_equal(arr,0)])
    arr = -1 * log10(arr)

    f = plt.figure(3)
    ax = f.add_subplot('111')

    clines = open(cfg.dataPath('daniel/heatplots/net{0}_top1676_communities.txt'.\
                                    format(net))).readlines()
    n_per_modules = [len(c.split('\t')) for c in clines]
    glines = open(cfg.dataPath('daniel/heatplots/net{0}_goterm_counts.txt'.\
                                    format(net))).readlines()
    n_per_go = dict([c.split('\t') for c in glines if c.strip() != ''])
    for k, v in n_per_go.iteritems(): n_per_go[k] = int(v.strip())

    big_mods = nonzero(greater(n_per_modules,min_module_size))[0]
    big_gos = set([ k for k, v in n_per_go.iteritems() if v > min_go_size])
    
    col_tits = [s.strip()  for s in l0.split('\t')[1:] if s.strip() != '']

    #FOR SOME WEIRD REASON, ONE OF THE COLUMNS THAT SHOULD BE A GO NAME IS 
    #CALLED 'V3'. AS V3 IS NOT PRESENT IN THE GO DESCRIPTIONS LIST,
    #I LEAVE IT OUT.

    acols = [ idx for idx, elt in enumerate(col_tits) if elt in big_gos ]
    arows = [ idx for idx, elt in enumerate(ids) if elt in big_mods ]

    arr = arr[arows][:, acols]

    thr = -1 * log10(array(prb_threshold))
    arr[greater(arr,thr[0])] = thr[0]
    arr[less(arr,thr[1])] = thr[1]
    
    arr = arr -  np.min(arr)
    arr = arr /  np.max(arr)

    go_modules = sum(arr, 0)
    final_cols = nonzero(less(go_modules,max_go_modules))[0]
    arr = arr[:,final_cols]

    return arr, \
        [ col_tits[idx] for idx in acols], \
        array([ ids[idx] for idx in arows])
