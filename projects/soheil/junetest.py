'''
Just a little routine to view the correlations of transcription factors and various genes. 
Written to examine the quality of network TF-TG edges.

Will probably never be used again...
'''
import matplotlib.pyplot as plt
import numpy as np
import cb.utils.colors as mycolors
from numpy import *
import scipy.io as sio
import cb.config as cfg
import cb.p.network.io as nio
names = ['FBgn0000411','FBgn0001325','FBgn0001180','FBgn0003900']


def run(num = 2):
    dfile = sio.loadmat(cfg.dataPath('soheil/expression_c4d_n4_intercluster.mat'))


    
    trgs, tfs = nio.getNet()
    bdgenes = nio.getBDTNP()
    bdset = set(bdgenes.keys())


    xs, ys, colors, corrs,lcorrs = [[] for i in range(5)]
    count = 0
    for k, v in bdgenes.iteritems():
        count += 1
        if count < num:
            continue
        if not trgs.has_key(k): continue
        trg = trgs[k]
        fsub = set(tfs.keys()).intersection(bdset)

        gexpr = bdgenes[k]['vals'][::50,4].flatten() #squeeze(dfile[k]) 
        fexpr = [bdgenes[fname]['vals'][::50,4].flatten() for fname in fsub]#[squeeze(dfile[fname]) for fname in fsub]
        
        
        print shape(fexpr)
        if len(fexpr )< 3: continue
        ct = mycolors.getct(len(fexpr))
        for idx, f in enumerate(fexpr):
            c = corrcoef(f, gexpr)[0,1]
            if not isfinite(c): c = 0
            lc = corrcoef(log(f), log(gexpr))[0,1]
            if not isfinite(lc): lc = 0
            corrs.append(c)
            lcorrs.append(lc)
            ys.append(gexpr)
            xs.append(f)
            colors.append([ct[idx]]* len(f))
        break
        if len(xs) > 10000:
            break
    
    cbest = argsort(-1 * abs(array(corrs)))
    

    f = plt.figure(1)
    f.clear()
    ax = f.add_subplot(111)
    inds = argsort(gexpr)

    for idx in cbest[:3]:
        import scipy.signal as ss
        import cb.utils.sigsmooth as sgs
        #k = sgs.gauss_kern(15)[8,:].flatten()
        #xconv = ss.convolve(xs[idx][inds],k) 
        
        xv = ss.medfilt(xs[idx][inds],1)
        yv = ys[idx][inds]
        print corrcoef(xv,yv)[0,1]
        print corrs[idx]
        
        ax.plot(ss.medfilt(xs[idx][inds],1), linewidth = 10, color = colors[idx][0])
        ax.plot(ys[idx][inds], linewidth = 10)

    #    print corrs[idx]
    #    ax.scatter(log(xs[idx]),log(ys[idx]),7, color =colors[idx], alpha = .3)
