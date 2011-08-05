#!/usr/bin/env python
from numpy import *
import numpy as np
import scipy.io as sio
import matplotlib as mpl
na = newaxis

import cb.p.network.io as nio
import cb.p.network.bdtnp.exp as exp
import cb.p.network.bdtnp.parser as bdparse

import cb.utils.plots as myplots
import cb.utils.colors as mycolors
import cb.utils.memo as mem
import cb.config as cfg


import textwrap as tw, itertools as it

def fetch_cluster_results(keys):
    cells = []
    for k in keys:
        try:
            res = sio.loadmat(cfg.dataPath('soheil/ben_results/cluster_results/results/'+\
                                               'expression_c4d_n4_t{0}.mat'.format(k)))
            
            cells.append((k, [  dict(gene = int(squeeze(r[0][0])),
                                     sign = r[1][0],
                                     module = tuple(sorted(r[2][0])),
                                     circuit = r[3][0],
                                     score = r[4][0])
                                for r in res['network_cell']]))
        except Exception, e:
            cells.append((k,[]))
            
    return cells

def cluster_elps(coords, tis, switch_axes = False):
        
    tis_means = array([np.mean([coords[:,idx[0]] for idx in g],0)  
                 for k, g in it.groupby(enumerate(tis),
                                        key = lambda x: x[1])])
    tis_vars = [ np.var([coords[:,idx[0]] for idx in g],0)
                 for k, g in it.groupby(enumerate(tis),
                                        key = lambda x: x[1])]
    centered =array( [coords[:,i] - tis_means[tis[i],:] for i in range(len(tis))])

    tis_cov = [sum([ dot(centered[idx[0],:][:,newaxis],
                         centered[idx[0],:][newaxis,:]) for idx in g],0) 
               for k, g in it.groupby(enumerate(tis),
                                      key = lambda x: x[1])]   
    tis_eigs = [linalg.eig(x) for x in tis_cov]

    import matplotlib.patches as patch

    all_elpses = []
    for j in range(2):
        elpses = []
        for i, elt  in enumerate(zip(tis_means, tis_vars, tis_eigs)):
            m, v ,e= elt
            vals, vecs = e
            
            if switch_axes:
                mean = m[:2] if j == 0 else [ m[0],m[2]]
                vars = array(vals[:2]) if j == 0 else array([vals[0],vals[2]]) 
                vecs  =  array(vecs[:2]) if j == 0 else array([vecs[0],vecs[2]]) 
            else:
                mean = m[:2] 
                vars = array(vals[:2])
                vecs  =  array(vecs[:2])
                
                
            vars = sqrt(vars)
            el = patch.Ellipse(mean,*( vars),\
                                   angle = 1.*180  / pi *\
                                   arctan(vecs[0][1]/vecs[0][0])) 
            elpses.append(el)
        all_elpses.append(elpses)
    return all_elpses


def get_coords( axes = 'gene', 
                rows = None,
                time_val = None,
                spatial_idxs = None,
                ids = None):
    bdnet = nio.getBDTNP()
    gene_matrix = array([v['vals'][:,time_val] 
                         for v in bdnet.values() 
                         if str(time_val + 1) in v['steps']])
    gene_matrix_keys = [k 
                for k in bdnet.keys() if str(time_val +1) in v['steps']]

    if axes == 'gene':
        import scipy.sparse as ssp
        import scipy.sparse.linalg as las
        import scipy.sparse.lil as ll
        adj = ssp.csr_matrix(gene_matrix.T)
        n_c = 3
        U,s, Vh = svd = las.svd(adj, n_c)
        filtered_genes = ll.lil_matrix(U)*ll.lil_matrix(diag(s)) *ll.lil_matrix(Vh)        
        xs_gene  = U[ids,0]
        ys_gene  = U[ids,1]
        zs_gene  = U[ids,2]
        
    elif axes == 'space':
        space_space =array([[ [r[idxs]  for idxs in sidxs]
                              for sidxs in spatial_idxs] 
                            for r in rows])
        space_space  = space_space[:, : , time_val]
    
        xs_gene = space_space[ids, 0]
        ys_gene = space_space[ids, 1]
        zs_gene = space_space[ids, 2]
    return xs_gene, ys_gene, zs_gene

def get_results(**kwargs):
    def set_results(**kwargs):
        cells = fetch_cluster_results([t[0] for t in kwargs.get('tsrt')])

        mod_list = list( set([m['module'] 
                          for c in cells
                          for m in c[1]] ) )
        mods =dict([(mod, [{'tissue':c[0],
                            'gene':m['gene']}
                           for c in cells 
                           for m in c[1] 
                           if m['module'] == mod] )
                    for mod in mod_list])
        tf_list = set(it.chain(*mods.keys()))
        gene_list = set([elt['gene'] for v in mods.values() for elt in v ])
        tfs = dict([(tf, [{'tissue':elt['tissue'],
                           'module':k,
                           'gene':elt['gene']}
                          for k, v in mods.iteritems()
                          if tf in k
                          for elt in v ])
                    for tf in tf_list])
        genes = dict([(g, [{'tissue':elt['tissue'],
                        'module':k}
                       for k, v in mods.iteritems()
                       for elt in v 
                       if elt ['gene'] == g])
                  for g in gene_list])
        return mods, genes, tfs
    return mem.getOrSet(set_results, 
                        **mem.rc(kwargs,
                                 on_fail = 'compute'))

def get_clusters(cluster_id = 4, 
                 switch_axes = False):
    all_members, ct_data = exp.recall_c2()
    all_members = array(all_members)
    c = all_members[cluster_id]
    c_unq = set(list(c))
    tissues = dict([('t_{0}'.format(i) , 
                     dict(cts = ct_data[equal(c,elt)]))
                    for i, elt in enumerate(c_unq)])
    bdnet = nio.getBDTNP()
    gene_cols, misc_cols, rows, row_nns = bdparse.read()
    spatial_idxs = array([misc_cols['x']['idxs'],
                     misc_cols['y']['idxs'],
                     misc_cols['z']['idxs']])

    #CHOOSE A TIME 
    time_val = cluster_id
    
    #ORGANIZE A BUNCH OF DATA.
    tsrt = sorted([(k,v) for k, v in  tissues.iteritems()],
                  key = lambda x: len(x[1]['cts']))[::-1]
    mods, genes, tfs = get_results(tsrt = tsrt, 
                                   name = 'cluster_{0}'.format(cluster_id))
    mod_srt = sorted(mods.iteritems(), key =lambda x: len(x[1]))[::-1]
    mod_ofinterest = mod_srt[0]


    f = myplots.fignum(3, (14,8))    
    ids,tis, cs = [[] for i in range(3)]
    n_tis = 60; ct = mycolors.getct(n_tis)
    for i,t in enumerate(tsrt[:n_tis]):
        nucs = [x[0] for x in t[1]['cts'] if x[1] == time_val][::5]
        ids.extend(nucs)
        tis.extend([i] * len(nucs))
        cs.extend([ct[i] for j in range( len(nucs))])
        
    #GET COORDS
    axes = 'gene'
    coords = array(list(get_coords(axes = axes,
                                   time_val = time_val,
                                   spatial_idxs = spatial_idxs,
                                   rows = rows,
                                   ids = ids)))
    #GET ELPSES
    all_elps = cluster_elps(coords, tis)
    module_scores = {}
    mods_of_interest = mod_srt[:3]
    for m in mods_of_interest:
        mcounts = []
        for i, e in enumerate(all_elps[0]):
            t = tsrt[i]
            tkey = t[0]
            mcounts.append( len([e 
                                 for elt in m[1] 
                                 if elt['tissue'] == tkey]))
        mvals = array(mcounts, float) / max(mcounts)
        module_scores[m[0]] = mvals

        
    #PLOT EM
    ax0 = f.add_subplot(121)
    ax1 = f.add_subplot(122)
    all_axes = [f.add_subplot(ax_str(ax0,ax1,ax2,ax3]

    cnv = mpl.colors.ColorConverter()
    for j, ax in enumerate([ ax0, ax1]):
        if switch_axes: myplots.padded_limits(ax, coords[0,:], coords[j+1,:])
        else:  myplots.padded_limits(ax, coords[0,:], coords[1,:])
        elps = all_elps[j]
        for i, e in enumerate(elps):
    
            
            yellow =squeeze(\
                    mpl.colors.rgb_to_hsv( array(cnv.to_rgb('yellow'))[na, na, :]))
            red = array([0.,1,1])
            mscore = module_scores.values()[j][i]
            hsv = yellow * (1 - mscore) + red * mscore
    
            e.set_alpha(.75)
            e.set_zorder(mscore)
            color = squeeze(mpl.colors.hsv_to_rgb(hsv[na,na,:]))
            e.set_facecolor(color)
            e.set_edgecolor('black')
            ax.add_patch(e)
            

    for annote_axis , plane in zip(*[[ax0, ax1],['XY','XZ']]):
        annote_axis.set_title('\n'.join(tw.wrap('''PCA projection in gene space of'''+\
                                        ' blastoderm nuclei for time = {0}.'.format(time_val)+\
                                        'Colors represent clusters used in'+\
                                        'model building\n{1} Axes, {2} Plane'.\
                              format(time_val,axes, plane),50)))
        annote_axis.set_xticks([])
        annote_axis.set_yticks([])    
    f.savefig(myplots.figpath('first_filtering_time{0}_Axis={1}.pdf'.\
                                  format(time_val,axes), delete = True))

    raise Exception()

def mod_list():
    path = cfg.dataPath("soheil/results_new/to ben-july 21/modules/denovo.mat")
    mods_d = sio.loadmat(path)['denovo']
    
    gene_mods = {}
    tf_mods = {}
    for m in mods_d:
        m = [tuple(list(elt[0])) for elt in m]
        gm = gene_mods.get(m[1], [])
        tm = tf_mods.get(m[0], [])
        gm.append(m[0])
        tm.append(m[1])
        gene_mods[m[1]] = gm
        tf_mods[m[0]] = tm

    


    return(tf_mods, gene_mods)
