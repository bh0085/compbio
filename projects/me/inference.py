import parse as wp
from numpy import *
import numpy as np
import itertools as it
import cb.utils.plots as myplots

import utils as wu
import cb.utils.memo as mem

def get_easy0(**kwargs):
    '''an easy inference using the arbitrary distance cutoff of 3000 bases
    and grabbing the n highest scoring edges globally'''
    def set_easy0(**kwargs):
        atype = kwargs.get('atype')
        simple = wp.get_simple_thr(atype = atype,
                                   dthr = 1500,
                                   dsign = -1,
                                   sthr = 1e-4)
        raise Exception()
        score_soft_cut = -136
        score_hard_cut = -90
        sids = wu.symbol_ids()
        prop_tuples = []
        for tf ,props in simple.iteritems():
           #for now, remove tfs that are not mappable
           if not tf in sids.keys(): continue
           ssrt = argsort(props['scores'])
           lscores = log10(props['scores'][ssrt])
       
           easy = nonzero(less(lscores , score_soft_cut))[0]
           medium = nonzero(greater(lscores, score_soft_cut)*\
                                less(lscores, score_hard_cut))[0]
           #generous_edges = concatenate([easy,medium])
           prop_tuples.append([(tf, props['genes'][g] ,
                                -(score_hard_cut - lscores[g]) \
                                    /(score_soft_cut - score_hard_cut)) 
                               for g in medium])
           prop_tuples.append([(tf, props['genes'][g],1) for g in easy])
           
        edgelist = array(list(it.chain(*prop_tuples)))
        edges = [( sids[e[0]] ,e[1].qualifiers['db_xref'][1][9:],e[2] )
                 for e in edgelist]
        return edges

    return mem.getOrSet(set_easy0, **mem.rc(kwargs))

def get_expression(**kwargs):
    atype = 'wormtile'
    simple = wp.get_simple_thr(**mem.rc(kwargs,
                                        atype = atype,
                                        dthr = 3000,
                                        dsign = -1,
                                        sthr = 1e-5))
    
    idfun = lambda g: g.qualifiers['db_xref'][1][9:] 
    genes = dict([(k,set([g
                          for g in v['gnames'] ])) 
                  for k,v in simple.iteritems()])
    gene_union = set([])
    for glist in genes.values(): gene_union.update(glist)

    
    crofs = wp.chromosome_offsets()
    gene_info = wp.gene_info(**mem.rc(kwargs))

    genome_coords = dict([(k,gene_info[k]['genomestart']) for i, k in enumerate(gene_union)])
    gene_idxs = dict([(k,i) for i, k  in enumerate(gene_union)])
    gnames = list(gene_union)
    gene_srtlist = argsort([e[1]  for e in sorted(gene_idxs.iteritems(), 
                                                  key = lambda x: genome_coords[x[0]])])
    gene_srtidxs= dict([(k,gene_srtlist[e]) for k ,e in gene_idxs.iteritems()]) 
    

    na = len(simple)
    ng = len(gene_union)
    gc = zeros((na,ng))

    
    assay_coords = dict([(k, i) for i , k in enumerate(simple.keys())])
    gene_counts =dict( [ (k,dict([(k2,len(list(g2)))  
                          for k2,g2 in it.groupby(
                            sorted(v['gnames']))])) 
                         for k,v in simple.iteritems()])
    for k,v in gene_counts.iteritems():
        for k2, v2 in v.iteritems():
            gc[assay_coords[k], gene_srtidxs[k2] ] = v2
    return 1*greater(gc,0)
def get_coexpression(gc, **kwargs):
    gc = gc/sum(gc, 0)
    f = myplots.fignum(3, (8,8))
    ax = f.add_subplot(111)
    
    for c in gc[:1]:
        #cplot = nonzero(greater(c,-1))[0][::10]
        #xs = array([genome_coords[k[0]]
        #            for k in sorted(gene_srtidxs.iteritems(),
        #                            key = lambda x: x[1])])
        ys = c

        #ax.scatter(*array([[gene_srtidxs[k],genome_coords[k]] for k in gene_union]).T)
        #ax.plot(xs, ys + random.rand(len(ys))*.1, color = random.rand(3),
        #        alpha = .25)

        
    cc = corrcoef(gc.T)
    ax.imshow(cc[:1000:1,:1000:1], aspect = 'auto')

    f.savefig(myplots.figpath('coex_counts_per_tissue.pdf'))
    
    return cc
    return genes,gene_counts, gene_info

def process_ccs(cc_in):
    cc = cc_in[:300,:300]
    for i,c in enumerate(cc): c[i] = 0
    import rpy2.robjects as robjects
    from rpy2.rlike.container import TaggedList
    from rpy2.robjects.packages import importr
    r = robjects.r

    base = importr('base')
    # create a numerical matrix of size 100x10 filled with NAs
    nc = nr = shape(cc)[0]

    from rpy2.robjects.numpy2ri import numpy2ri


    m = numpy2ri(cc)#robjects.r['matrix'](v, nrow = nr, ncol = nc)
    
    

    biclust = importr("biclust")
    mb = biclust.binarize(m,.90)

    #hcv = r.hclust(r.dist(mb))
    #hcv = r.hclust(r.dist(mb))
    #hm = r.heatmap(mb)

    #raise Exception()

    out = biclust.biclust(m, method=biclust.BCPlaid())
    n_bc = out.do_slot('Number')
    rows = array(out.do_slot('RowxNumber'))
    cols = array(out.do_slot('NumberxCol')).T

    return rows, cols, array(m), array(mb)

def pvclust(exprs):
    
    pass


def process_rc(cc,rows,cols, meth = 'binary'):
    rmembers = zeros(len(rows)) + 4
    cmembers = zeros(len(rows)) + 4
    for i,r in enumerate(rows):
        rmembers[i] = argmax(r) if np.max(r) > 0 else len(r)
    for i,c in enumerate(cols):
        cmembers[i] = argmax(c) if np.max(c) > 0 else len(c)
    
    rorder = argsort(rmembers)
    corder = argsort(cmembers)
    
    f = myplots.fignum(3,(8,8))
    ax = f.add_subplot(111)
    ax.imshow(cc[rorder][:,corder], aspect = 'auto', interpolation = 'nearest')
    
    f.savefig(myplots.figpath('biclustered_expr_{0}.pdf'.format(meth)))
    raise Exception()

