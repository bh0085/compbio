import scipy.io as sio
from numpy import *
import compbio.config as cfg
import os
import matplotlib.pyplot as plt
import compbio.utils.bsub_utils as butils
import compbio.utils.plots as myplots
def view():

    fnums = range(3455)[::]
    #outputs =[ load_data( open(cfg.dataPath('batch/tmp/run_mcmc_{0:05}_tmp001.mat'.\
    #  
    outputs =[ sio.loadmat(cfg.dataPath('batch/tmp/mcmc_{0:05}_tmp001.mat'.format(num)))  for num in fnums ]

    douts = []
    for output in outputs:
        try:
            o00 = output['out_struct'][0][0]
            dout = dict([(k, o00[i]) for i , k in enumerate([elt[0] for elt in o00.dtype.descr])])
            douts.append(dout)
        except Exception, e:
            continue
    
    ss, ir = array([(squeeze(o['stay_same']),squeeze(o['improve_ratio'])) for o in douts]).T 
    ss += random.rand(*shape(ss))/100
    ir += random.rand(*shape(ir))/100

    f = myplots.fignum(1, (8,8))
    f.clear()
    ax = f.add_subplot(111)
    ax.set_xlabel('Stay Same')
    ax.set_ylabel('Improve Ratio')
    plt.scatter(ss, ir, 5)
    
    
    
def view2():
    files = [l for l in os.listdir(cfg.dataPath('batch/outputs')) if 'mcmc' in l]
    ids = [l[0:10] for l in files]
    ids = ids[::10]

    inps = [ butils.load_data(i,'input') for i in ids ]
    outs = [ butils.load_data(i,'output') for i in ids ]
    
    idxs_good = nonzero(greater([elt.get('improve_ratio') for elt in outs], .2) *
                        greater([elt.get('beta') for elt in inps], 4) )[0]
    
    outs = [ o for  idx, o in enumerate(outs) if idx in idxs_good]
    inps = [ i for  idx, i in enumerate(inps) if idx in idxs_good]

    params = inps[0].keys()
   
    f =myplots.fignum(1, (8,8))
    
    params = params


    for i, p in enumerate(params):
        ax = f.add_axes([.05, i * ( 1./len(params)) ,  .9, 1./len(params)] , title = p)
        #ax.set_yticks([])
        #ax.set_xticks([])

        xvals = [elt.get(p) for elt in inps]
        if type(xvals[0]) == str: continue
        yvals = [elt.get('improve_ratio') for elt in outs] 
        yvals2 = [elt.get('stay_same') for elt in outs] 
        
    

        yvals += random.rand(*shape(yvals)) * (max(yvals) - min(yvals))/50
        yvals2 += random.rand(*shape(yvals)) * (max(yvals) - min(yvals))/50
        xvals += random.rand(*shape(xvals)) * (max(xvals) - min(xvals))/50
        ax.scatter(xvals,yvals)
        
        #ax.scatter(xvals , yvals + yvals2,   25, color = 'red')
        ax.annotate(p, [0,0], xycoords = 'axes fraction', ha = 'left', va = 'bottom')
    
    f.savefig(cfg.dataPath('figs/soheil/broad_run0_psplits.ps'))

    return inps
    
