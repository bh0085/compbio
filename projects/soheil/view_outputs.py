import scipy.io as sio
from numpy import *
import compbio.config as cfg
import os, re
import matplotlib.pyplot as plt

import compbio.utils.bsub_utils as butils
import compbio.utils.plots as myplots
import cb.utils.seismic as seismic
import itertools as it
import compbio.utils.memo as mem

figtemplate = cfg.dataPath('figs/soheil/vo_{0}.pdf')

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
    

def fetch_genes():
    files = [l for l in os.listdir(cfg.dataPath('batch/tmp')) if 'mcmc' in l]
    fpaths = [os.path.join(cfg.dataPath('batch/tmp') ,f) for f in files]
    ids = [l[0:10] for l in files]


    inps = [ butils.load_data(i,'input') for i in ids ]
    l_info= {}
    for l, elt in enumerate(zip(fpaths, inps)):
        f, inp = elt
        if inp['out_iter_num'] == 2: continue
        print inp['filename']
        clustname = re.search(re.compile('_([^_]+)\.mat'),inp['filename']).group(1)
        l_info[l] = {}
        l_info[l]['cname']= clustname
        l_info[l]['filename'] = inp['filename']

        data = sio.loadmat(f)
        l_info[l]['stay_same'] = data['stay_same']
        l_info[l]['improve_ratio'] = data['improve_ratio']
        l_info[l]['error_test'] = data['error_test']
        l_info[l]['error_test'] = data['error_test']
    



def view4():
    
    files = [l for l in os.listdir(cfg.dataPath('batch/tmp')) if 'mcmc' in l]
    fpaths = [os.path.join(cfg.dataPath('batch/tmp') ,f) for f in files]
    ids = [l[0:10] for l in files]
    inps = [ butils.load_data(i,'input') for i in ids ]

    idxs_good = nonzero(greater([elt.get('out_iter_num') for elt in inps], -1 ))[0]
    inps = [inps[i] for i in idxs_good]
    fpaths = [fpaths[i] for i in idxs_good]

    
    termgroups, cnames, xvals, gvals, yvals, colors, rads,tfs, all_coefs  =\
        [],[],[], [], [], [], [], [], []
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
        l_info[l]['error_test'] = data['error_test']

        import compbio.utils.colors as mycolors
        ct = mycolors.getct(len(data['gene_names']))
        
        term_list = [ list(it.chain(*mod)) for mod in data['model']]
        fac_list = [list(it.chain(*t)) for t in term_list]

        seen = set()
        all_coefs.append(data['coefs_dic_nonlinear'])
        coefs = data['coefs_dic_nonlinear']
        nlcof_all = open(cfg.dataPath(\
                'network/network_predmodel/regressionwts/nonlinear_all/nw_{0}.sif'.\
                    format(l)),'w')
            
        nlcof_sing = open(cfg.dataPath(\
                'network/network_predmodel/regressionwts/nonlinear_sing/nw_{0}.sif'.\
                    format(l)),'w')

        tfnames = data['tf_names']
        tgnames = data['gene_names']
        

        for i, terms in enumerate(term_list):
            if i in (5,49,53,30,17,8,38) :
                if sum(terms) > 0:
                    raise Exception()
            terms = [t -1 for t in terms]
            for j, term in enumerate(terms):
                if len(term) == 1:
                    wt = coefs[i][0][0][j]
                    nlcof_sing.write('{0}\t{1}\t{2}\n'.\
                                         format(tfnames[term][0][0],
                                                tgnames[i][0], wt))
                
                
                for k, fact in enumerate(list(set(term))):
                    wt = coefs[i][0][0][j]
                    nlcof_all.write('{0}\t{1}\t{2}\n'.\
                                     format(tfnames[fact][0][0],
                                            tgnames[i][0][0],
                                            wt))
                    
                    gvals.append([i] * (len(term)+1))
                    yvals.append([fact] *( len(term)+1))
                    colors.append([ct[c] for c in sorted(term)]+[1,1,1])
                    tfs.append([c for c in sorted(term)])
                    rads.append(((arange(1,len(term) +2 )**2) * 50)[::-1])
                    xvals.append([l] * (len(term)+1))
                
        nlcof_all.close()
        nlcof_sing.close()
    
    return cnames, xvals, gvals, yvals, colors, rads, l_info, tfs, coefs

def score_genes(cnames, xvals, gvals, yvals, colors, rads,l_info):
    errors = array([mean(x['error_test'],1) for x in l_info.values()])
    errmean = mean(errors,0)
    g_errsrt = argsort(errmean)
    g_errsrt = [g for g in g_errsrt if errmean[g] != 0]
    return g_errsrt

def errors():
    files = [l for l in os.listdir(cfg.dataPath('batch/tmp')) if 'mcmc' in l]
    fpaths = [os.path.join(cfg.dataPath('batch/tmp') ,f) for f in files]
    ids = [l[0:10] for l in files]


    inps = [ butils.load_data(i,'input') for i in ids ]

    idxs_good = nonzero(greater([elt.get('out_iter_num') for elt in inps], -1 ))[0]
    inps = [inps[i] for i in idxs_good]
    fpaths = [fpaths[i] for i in idxs_good]

    errors, staysames, improves = [], [], []
    for l, elt in enumerate(zip(fpaths, inps)):
        f, inp = elt

        data = sio.loadmat(f)
        errors.append(data['error'])
        staysames.append(data['stay_same'])
        improves.append(data['improve_ratio'])
        gnames = data['gene_names']
        
    return errors, staysames, improves, gnames
def names():

    files = [l for l in os.listdir(cfg.dataPath('batch/tmp')) if 'mcmc' in l]
    fpaths = [os.path.join(cfg.dataPath('batch/tmp') ,f) for f in files]
    data =  sio.loadmat(fpaths[0])

    gnames = data['gene_names']
    tfnames = data['tf_names']
    return gnames, tfnames

def show_errors( errors, staysames, improves, gnames):
    figtitle = 'show_errors'
    f = myplots.fignum(3, (12,6))
    ax = f.add_axes([.05,.05,.25,.9])

    import scipy.signal as ss
    for all_errs in errors[0:1]: 
        for e in all_errs.flatten()[:5]:
            ax.plot(ss.medfilt(e.flatten()**2,51))

    get_worse = 1 - (array(staysames) + \
                         array(improves))
    ax2 = f.add_axes([.3,.05,.65,.9])
    seismic.seismic(squeeze([get_worse, staysames, improves]), stacked = True,
                    colors = [[1,0,0], [0,0,0], [0,0,1]],
                    ax = ax2, linewidth = 10, label_y = False)
    
    f.savefig(figtemplate.format(figtitle))

def view4_show0(cnames, xvals, gvals, yvals, colors, rads, l_info, gnum = 59):
    
    seen = set()
    offs_mag = .3
    xofs,xnew,yofs, ynew = [], [], [], [] #[xv for xv in xvals], [yv for yv in yvals]
    for v in zip(xvals, yvals, gvals):

        xy0 = v[0][0], v[1][0]
        xy = tuple([x for x in xy0])
        #check to see if the current xy has been seen.
        #if so increment until unique.
        while 1:
            if xy in seen:
                xy = tuple([elt + offs_mag for elt in xy])
            else:
                break
        
        xnew.append([xy[0] for i in range(len(v[0]))])
        ynew.append([xy[1] for i in range(len(v[0]))])

        xofs.append([xy[0] - xy0[0] for i in range(len(v[0]))])
        yofs.append([xy[1] - xy0[1] for i in range(len(v[0]))])

        if v[2][0] != gnum:
            continue
        else:
            seen.add(xy)




    g_equal = nonzero(equal([x for x in it.chain(*gvals)],gnum))[0]
    if len(g_equal) == 0:
        print 'G {0} appears to not be in the list'.format(gnum)
    gset = set(g_equal)

    xvals_old = xvals
    yvals_old = yvals
    xvals = xnew
    yvals = ynew
    
    xvals = array(list(it.chain(*xvals)))[g_equal]
    yvals = array(list(it.chain(*yvals)))[g_equal]
    
    xvals_old = array(list(it.chain(*xvals_old)))[g_equal]
    yvals_old = array(list(it.chain(*yvals_old)))[g_equal]
    
    xofs = array(list(it.chain(*xofs)))[g_equal]
    yofs = array(list(it.chain(*yofs)))[g_equal]

    colors = array(list(it.chain(*colors)))[g_equal]
    rads = array(list(it.chain(*rads)))[g_equal]
    
    


    vecs = zeros((max(xvals_old) + 1, max(yvals_old)+1))
    for x, y in zip(xvals_old, yvals_old):
            vecs[x,y] = 1.
        
    #import hcluster
    #clusters = hcluster.fclusterdata(vecs,.1,criterion='inconsistent',method = 'complete' )
    import mlpy
    HC = mlpy.HCluster(method='euclidean', link='complete')
    clusts = HC.compute(vecs )
    k = 15
    cut = HC.cut(HC.heights[-k]) 
    cut_s = sort(cut)
    crank = argsort(argsort(cut))

    fig = myplots.fignum(3, (35,15))
    ax =fig.add_axes([0,0,1,1])

    clst_mems = [cut_s[c] for c in crank]
    clust_colors = ones(3) * linspace(.2,.9,len(set(cut)))[:, newaxis]
    
    ax.scatter(crank, ones(len(crank) )* .5, 1000, color = clust_colors[clst_mems])
    
    ax.scatter(crank[xvals_old] +xofs, 
               yvals,rads, 
               color = colors)


    ax.annotate('Functional motifs for gene: {0}\nIn {1} clusters'.\
                    format('gene',100),
                [0,1], xycoords = 'axes fraction',va = 'top'
                )

    
    figname = 'gene_{0}_motif_recurrence_circles'.format(gnum)
    fig.savefig(figtemplate.format(figname))


def find_clusters(tfs,cnames, xvals, gvals, yvals, colors, rads, l_info,
                  gnames, tfnames):
    

    

    gs = [x[0] for x in gvals]
    grps = [(k, list(g)) for k, g in it.groupby(\
            sorted(enumerate(gs),key = lambda x: x[1]), 
            key = lambda x:x[1]
            )]
    print '{0} genes'.format(len(grps))
    mlens, modules = {}, {}
    for kg,gg in grps:
        g_ids = [g[0] for g in gg]
        clusters = [(i,xvals[i][0]) for i in g_ids]
        cgrps = [(k,list(g)) for k, g in it.groupby(\
                sorted(clusters, key = lambda x: x[1]),
                key = lambda x:x[1])]
        tsets = []
        for k,cg in cgrps:
            c_ids = [c[0] for c in cg]
            terms = [tfs[i] for i in c_ids]
            tset = list(set([tuple(s) for s in terms]))
            tsets.append(tset)
    
        modules[kg] = [(k,list(g)) for k, g in it.groupby(\
                sorted([(cid,tuple(set(term)))
                        for cid, termlist in enumerate(tsets)
                        for term in termlist],
                       key = lambda x: x[1]),
                key = lambda x: x[1])]


    
    mgrps = dict([(mk,dict([(k, list(g)) for k, g in it.groupby(
            sorted(it.chain(*[[m_subelt[1] for m_subelt in m_elt[1]]
                              for m_elt in mv]))
            )])) 
            for mk,mv in modules.iteritems()])

    munq = dict([(mk,[mgelt for mgelt in mg.keys()]) 
                 for mk,mg in mgrps.iteritems()])
    


    mods = {}
    mod_genes = {}
    for mk, ms in munq.iteritems():
        for m in ms:
            mods[m] = mods.get(m, 0) + 1
            mod_genes[m] = mod_genes.get(m, [])
            mod_genes[m].append([gnames[mk], len(mgrps[mk][m])])
            
    

    mdubs = [ (k, v) for k, v in mods.iteritems() if len(k) == 2]
    mtrips = [ (k, v) for k, v in mods.iteritems() if len(k) == 3]
    
    trips_30 =sorted(mtrips, key = lambda  x: x[1] )[::-1][:30]
    gene_30 = [([tfnames[j][0][0] for j in i[0] ],
                [[g[0][0][0],g[1]] for g in mod_genes[i[0]]],
                i[1]) 
               for i in trips_30]    

    dubs_30 =sorted(mdubs, key = lambda  x: x[1] )[::-1][:30]
    dgene_30 = [([tfnames[j][0][0] for j in i[0] ],
                 [[g[0][0][0],g[1]] for g in mod_genes[i[0]]],
                 i[1]) 
               for i in dubs_30]
    #srt = sorted(it.chain(*munq))
    
    #mgrps = [[k, list(g)] for k, g in it.groupby(srt)]
    return gene_30, dgene_30
        
        
    raise Exception()


def setModules(**kwargs):
    files = [l for l in os.listdir(cfg.dataPath('batch/tmp')) if 'mcmc' in l]
    fpaths = [os.path.join(cfg.dataPath('batch/tmp') ,f) for f in files]
    ids = [l[0:10] for l in files]
    inps = [ butils.load_data(i,'input') for i in ids ]

    modules = {}
    lin_modules = {}
    for fidx, f in enumerate(fpaths):
        print 'Getting module info for: {0}'.format(f)
        data = sio.loadmat(f)
        tfnames = [d[0][0] for d in data['tf_names']]
        tgnames = [d[0][0] for d in data['gene_names']]
        coefs =   [d[0][0] for d in data['coefs_dic_nonlinear']]
        inp = inps[fidx]
        
        term_list = [ list(it.chain(*mod)) for mod in data['model']]
        for j,terms in enumerate(term_list):
            if sum([len(t) for t in terms]) == 0:continue
            for k, t in enumerate(terms):
                mod =   tuple([tfnames[i] for i in sorted(t - 1)])
                mod_d = modules.get(mod, dict(genes = [], 
                                              coefs = [],
                                              fpaths = [],
                                              clust_fpaths = []))
                mod_d['genes'].append(tgnames[j])
                mod_d['coefs'].append(coefs[j][k])
                mod_d['clust_fpaths'].append(inp['filename'])
                mod_d['fpaths'].append(f)
                modules[mod] = mod_d

        lin_coefs =   [d[0][0] for d in data['coefs_dic_nonlinear']]
        term_list = [ list(it.chain(*mod)) for mod in data['model_linear']]
        for j,terms in enumerate(term_list):
            if sum([len(t) for t in terms]) == 0:continue
            for k, t in enumerate(terms):
                mod =   tuple([tfnames[i] for i in sorted(t - 1)])
                mod_d = lin_modules.get(mod, dict(genes = [], 
                                                  coefs = [],
                                                  fpaths = [],
                                                  clust_fpaths = []))
                mod_d['genes'].append(tgnames[j])
                mod_d['coefs'].append(coefs[j][k])
                mod_d['fpaths'].append(f)
                mod_d['clust_fpaths'].append(inp['filename'])

                lin_modules[mod] = mod_d
    return modules,lin_modules
   
def modules(reset = False):
    return mem.getOrSet(setModules, 
                        **mem.rc({},
                                 reset = reset, 
                                 hardcopy = True,
                                 on_fail = 'compute'))


        
