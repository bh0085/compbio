import scipy.io as sio
from numpy import *
import compbio.config as cfg
import os, re
import matplotlib.pyplot as plt
import compbio.utils.bsub_utils as butils
import compbio.utils.plots as myplots
import itertools as it
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
    
    #idxs_good = nonzero(greater([elt.get('improve_ratio') for elt in outs],, .2 )[0]
    idxs_good = range(len(outs))
 
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
    raise Exception()

    return inps
    
def view3():
    
    files = [l for l in os.listdir(cfg.dataPath('batch/tmp')) if 'mcmc' in l]
    fpaths = [os.path.join(cfg.dataPath('batch/tmp') ,f) for f in files]
    ids = [l[0:10] for l in files]
    ids = ids[::10]

    inps = [ butils.load_data(i,'input') for i in ids ]
    idxs_good = nonzero(greater([elt.get('out_iter_num') for elt in inps],2 ))[0]
    inps = [inps[i] for i in idxs_good]
    fpaths = [fpaths[i] for i in idxs_good]


    fig = myplots.fignum(3, (35,15))
    ax =fig.add_axes([0,0,1,1])

    for f, inp in zip(fpaths, inps):
        if inp['out_iter_num'] == 2: continue
        print inp['filename']
        
        data = sio.loadmat(f)
        
        
        import compbio.utils.colors as mycolors
        ct = mycolors.getct(len(data['gene_names']))
        
        term_list = [ list(it.chain(*mod)) for mod in data['model']]
        fac_list = [list(it.chain(*t)) for t in term_list]

        xvals, yvals, colors, rads  = [], [], [], []
        for i, terms in enumerate(term_list):
            for j, term in enumerate(terms):
                for k, fact in enumerate(term):
                    xvals.extend([i] * len(term))
                    yvals.extend([fact] * len(term))
                    colors.extend([ct[c] for c in sorted(term)])
                    rads.extend(((arange(1,len(term) +1 )**2) * 50)[::-1])
                    
                   
        vecs = zeros((len(fac_list), len(fac_list)))
        for i, fl in enumerate(fac_list):
            for f in fl:
                vecs[i,f] = 1
        
        #plt.imshow(vecs)

        #ax1 = fig.add_subplot(121)
        #ax2 = fig.add_subplot(122)
        import hcluster
        clusters = hcluster.fclusterdata(vecs,1.1,criterion='inconsistent',method = 'complete' )
        
        #ax1.imshow(vecs)
        #ax2.imshow(vecs[argsort(clusters)])

        #raise Exception()
            
                    
        csrt = argsort(argsort(clusters))
        xvals2 = [ csrt[x] for x in xvals]

        #raise Exception()
        plt.scatter(xvals2, yvals,rads, color = colors)
        raise Exception()

    raise Exception()
    
def view4():
    
    files = [l for l in os.listdir(cfg.dataPath('batch/tmp')) if 'mcmc' in l]
    fpaths = [os.path.join(cfg.dataPath('batch/tmp') ,f) for f in files]
    ids = [l[0:10] for l in files]
    ids = ids[::10]


    inps = [ butils.load_data(i,'input') for i in ids ]

    idxs_good = nonzero(greater([elt.get('out_iter_num') for elt in inps], -1 ))[0]
    inps = [inps[i] for i in idxs_good]
    fpaths = [fpaths[i] for i in idxs_good]



    fig = myplots.fignum(3, (35,15))
    ax =fig.add_axes([0,0,1,1])

    
    cnames, xvals, gvals, yvals, colors, rads  = [], [], [], [], [], []
    l_info= {}
    for l, elt in enumerate(zip(fpaths, inps)):
        f, inp = elt
        if inp['out_iter_num'] == 2: continue
        print inp['filename']
        clustname = re.search(re.compile('_([^_]+)\.mat'),inp['filename']).group(1)
        cnames.append(clustname)
        l_info[l] = {}
        l_info[l]['cname']= clustname
        l_info[l]['filename'] = inp['filename']

        data = sio.loadmat(f)
        l_info[l]['stay_same'] = data['stay_same']
        l_info[l]['improve_ratio'] = data['improve_ratio']

        import compbio.utils.colors as mycolors
        ct = mycolors.getct(len(data['gene_names']))
        
        term_list = [ list(it.chain(*mod)) for mod in data['model']]
        fac_list = [list(it.chain(*t)) for t in term_list]

        for i, terms in enumerate(term_list):
            for j, term in enumerate(terms):
                for k, fact in enumerate(term):
                    gvals.extend([i] * len(term))
                    yvals.extend([fact] * len(term))
                    colors.extend([ct[c] for c in sorted(term)])
                    rads.extend(((arange(1,len(term) +1 )**2) * 50)[::-1])
                    xvals.extend([l] * len(term))
                    

        #plt.imshow(vecs)

        #ax1 = fig.add_subplot(121)
        #ax2 = fig.add_subplot(122)

        #ax1.imshow(vecs)
        #ax2.imshow(vecs[argsort(clusters)])

        #raise Exception()
       
    return cnames, xvals, gvals, yvals, colors, rads, l_info

    
def view4_show0(cnames, xvals, gvals, yvals, colors, rads, l_info)

    #csrt = argsort(argsort(clusters))
    #xvals2 = [ csrt[x] for x in xvals]
    gvals = array(gvals)
    gnum = 59
    g_equal = nonzero(equal(gvals,gnum))[0]
    xvals = array(xvals)[g_equal]
    yvals = array(yvals)[g_equal]
    colors = array(colors)[g_equal]
    rads = array(rads)[g_equal]
    
    
                   
    vecs = zeros((max(xvals) + 1, len(fac_list)))
    for x in xvals:
        for y in yvals:
            vecs[x,y] = 1
        
    import hcluster
    clusters = hcluster.fclusterdata(vecs,1.1,criterion='inconsistent',method = 'complete' )
                          
    csrt = argsort(argsort(clusters))
    xvals2 = [ csrt[x] for x in xvals]
 

    plt.scatter(xvals2 , yvals,rads, color = colors)
    ax.annotate('Functional motifs for gene: {0}\nIn {1} clusters'.\
                    format(data['gene_names'][gnum],100),
                [0,1], xycoords = 'axes fraction',va = 'top'
                )


    raise Exception()
