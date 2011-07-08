import matplotlib.pyplot as plt
from numpy import *
import numpy as np
import compbio.config as cfg
import compbio.utils.plots as myplots
import cb.utils.colors as mycolors

def load(net = 2, num = 1676,
         min_module_size = 10,
         min_go_size = 5,
         max_go_modules = 2,
         prb_threshold = [1e-6,.01]):
    fopen = open(cfg.dataPath('daniel/heatplots/net{0}_top{1}_heatplot_matrix.txt'.format(net,num)),'r')
    l0 = fopen.readline()

    arr = array([[float(elt) for elt in l.split('\t')] for  l in  fopen.xreadlines() if l.strip() != ''])
    ids = arr[:,:1]
    arr = arr[:,1:]

    arr[equal(arr,0)] = np.min(arr[not_equal(arr,0)])
    arr = -1 * log10(arr)


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

    acols = array([ idx for idx, elt in enumerate(col_tits) if elt in big_gos ])
    arows = array([ idx for idx, elt in enumerate(ids) if elt in big_mods ])

    arr = arr[arows][:, acols]

    thr = -1 * log10(array(prb_threshold))
    arr[greater(arr,thr[0])] = thr[0]
    arr[less(arr,thr[1])] = thr[1]
    
    arr = arr -  np.min(arr)
    arr = arr /  np.max(arr)

    go_modules = sum(arr, 0)
    final_cols = nonzero(less(go_modules,max_go_modules)*\
                             greater(go_modules,0))[0]
    acols = acols[final_cols]
    arr = arr[:,final_cols]

    return arr, \
        array([ col_tits[idx] for idx in acols]), \
        array([ ids[idx] for idx in arows])

def srt_heatmap(net = 3,
                all_module = False):
    import compbio.projects.bsort.bsort as bs

    arr, cols, rows = load(net = net,
                              max_go_modules = 15,
                              min_go_size = 5, 
                              min_module_size = 10)

    arr2_510, srts = bs.run0(arr = arr, itr = 1, meth = 'moment')
    arr2_510 = arr2_510[:,::-1]
    csrts = [s for s in srts if len(s) == len(cols)]
    rsrts = [s for s in srts if len(s) == len(rows)]

    c0 = array(cols)
    r0 = array(rows)
    for c in csrts: cols = cols[c]
    for r in rsrts: rows = rows[r]

    fopen =open( cfg.dataPath('daniel/heatmaps_sorted/hm_net{0}.txt'.format(net)), 'w')
    fopen.write('FORMAT: L1 :GO Terms (Columns), L2: Modules (Rows), L3+ Pvals thresholded between .01, .001\n')
    fopen.write('\t'.join([str(elt) for elt in cols]) + '\n')
    fopen.write('\t'.join([str(squeeze(elt)) for elt in rows]) + '\n')

    dmat = arr2_510
    for row in dmat:
        fopen.write('\t'.join(['{0}'.format(elt) for elt in row])+'\n')
    fopen.close()

    f = myplots.fignum(3, (8,40)) if net == 3 else myplots.fignum(3, (8,10))
    ax = f.add_axes([.4,.05,.55,.9], aspect = 'auto')
    
    goterms = [str(elt) for elt in cols]

    if not all_module:
        idx_omit = goterms.index('all')
        arr2_510 = vstack((arr2_510[:idx_omit],arr2_510[idx_omit+1:]))
        goterms = goterms[:idx_omit] + goterms[idx_omit+1:]
        
        

    fopen = open(cfg.dataPath('daniel/go_accession_name_map.txt'))
    gotexts = {}
    for l in fopen.xreadlines():
        k,v = l.split('\t')
        gotexts[k] = v
    row_labels = []
    for g in goterms:
        if g in gotexts.keys(): row_labels.append(gotexts[g].strip())
        else: row_labels.append(g.strip())
    


    ax.set_yticks(arange(len(goterms))+.25)
    ax.set_yticklabels(row_labels, size = '4')
    ax.set_xticks([])
    ax.set_xlabel('Modules')

    ax.set_xticks(arange(len(rows))+.25)
    ax.set_xticklabels([str(int(r[0])) for r in rows],
                       rotation= 90,
                       size = 'small')

    

    cm = mycolors.blackbody(flip = True)
    im =ax.imshow(arr2_510[:,:] * 4 + 2,
              cmap = plt.get_cmap('OrRd'),
              aspect = 'auto', 
              interpolation = 'nearest'
              )
    plt.colorbar(im)
    f.savefig(cfg.dataPath('daniel/heatmaps_sorted/hm_net{0}_{1}.pdf'.\
                               format(net, 'with_all' if all_module else 'no_all')))

    
    
