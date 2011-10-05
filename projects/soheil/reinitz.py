import compbio.utils.plots as myplots
import compbio.utils.colors as mycolors
import os, itertools as it,re
import compbio.utils.memo as mem
import compbio.config as cfg

import cb.p.network.utils as nutils
import cb.p.nfilt.utils as nfutils
import cb.p.network.io as nio

from numpy import *
import numpy as np

def id_map(**kwargs):
    def set_id_map(**kwargs):
     fname = cfg.dataPath('reinitz/28-7-2011-1-56-6-30-0/txt/byGenes')
     
     gsums = open(cfg.dataPath('flybase/gene_summaries.tsv'))
     gmap = open(cfg.dataPath('flybase/gene_map.tsv'))
     gassoc = open(cfg.dataPath('flybase/gene_association.fb'))
     
     gname_orig =  [ os.path.splitext(f)[0].lower() for f   in  os.listdir(fname) ] 
     gfiles =dict(  [ (gname_orig[i], os.path.join(fname,f)) for i, f in  enumerate(os.listdir(fname)) ] )
     gname_map = dict([( re.sub( re.compile('[^a-z]'),'',g), g) for g in gname_orig])
     gnames = gname_map.keys()
     
     glines = dict([(k.lower(),[]) for k in gnames])
     
     lines_kept = {}
     for i, g in enumerate(gassoc.xreadlines()):
         if g[0] == '!': continue
         g0 = g
         g = re.sub( re.compile('[^a-z]'),'', g.lower().split('\t')[9].strip())
         for k,v in glines.iteritems():
     
             if k == g: 
                 v.append((i,g))
                 lines_kept[i] = g0
         
     
     matches = glines
     ids = {}
     for k, v in matches.iteritems():
         names =  [ l[1] for l in v] 
         line_nums =  [ l[0] for l in v] 
         these_ids = [lines_kept[i].split('\t')[1].strip() for i in line_nums] 
         #just hacking here... for sloppy paired I use the first id...
         #alas...
         ids[k] = tuple(sorted(set(these_ids)))[0]
     
     return dict([ (idval, {'file': gfiles[gname_map[k]], 'name':gname_map[k]}) for k, idval in ids.iteritems()])
         #name_grps = dict([(gpkey, list(g)) for gpkey, g in  it.groupby(sorted(names))])
         #print k
         #print [ (gk, len(gv)) for gk, gv in name_grps.iteritems()] 
    return mem.getOrSet(set_id_map,**mem.rc(kwargs,on_fail = 'compute'))

def datafiles(**kwargs):
    def set_datafiles(**kwargs):
        out ={}
        idmap = id_map(**mem.sr(kwargs))
        for k,v in idmap.iteritems():
            out[k] = array([ [float(e) for e in re.compile('\s+').split(l.strip())] for l in open(v['file']).readlines() if l[0] in '0123456789'])
        return out
    return mem.getOrSet(set_datafiles, **mem.rc(kwargs,
                                                on_fail = 'compute'))


def check_bdtnp(**kwargs):
    ids = id_map(**mem.sr(kwargs))
    
    btgs = nio.getBDTNP()
    stfs,stgs = nio.getNet()
    
    sxs = dict([(i, {}) for i in ids.keys()])
    for i, elt  in enumerate(ids.iteritems()):
        gid, gene_val = elt
        if stfs.has_key(gid) : sxs[gid]['stf'] = True;
        if btgs.has_key(gid) : sxs[gid]['btg'] = True;
        if stgs.has_key(gid) : sxs[gid]['stg'] = True;

    return sxs

def make_kde(vals):                    
    from scipy.stats import gaussian_kde
    count = len(vals)
    dist = array(vals)
    density = gaussian_kde(dist)
    density.covariance_factor = lambda : .5
    density._compute_covariance()
    kde = (density)
    
    return count, kde


import cb.p.nfilt.comparator as comp

def get_reinitz_data(**kwargs):

    ofs = kwargs.get('ofs',0)
    do_plot_coords = kwargs.get('plot_coords',False)
    do_plot_vals = kwargs.get('plot_vals',False)

    idm= id_map()
    df = datafiles(**mem.rc(kwargs))

    #I'm not sure exactly how this dataset works but
    #each nuclei has a bunch of numbers that appear to be
    #monotonically increasing.
    #
    #I just take the first instance.
    nums = dict([(k,v[:,0]) for k, v in df.iteritems()])
    nuc_count = len(set(nums.values()[2]))
   
    values = dict([(k,v[nuc_count *ofs: nuc_count *(ofs + 1),-1]) for k, v in df.iteritems()])
    coords = dict([(k,v[nuc_count *ofs :nuc_count *(ofs + 1),1:3]) for k, v in df.iteritems()])

    #to check the basic consistency of the data, enable the plot routines.
    #I suppose that I could do this for all of the nuclei occurences...
    #right now, only the first is used.
    if do_plot_coords:
        f = myplots.fignum(1,(8,8))
        ax = f.add_subplot(111)
        ct = mycolors.getct(len(values))
        for i,k in enumerate(values.keys()):
            ax.scatter(coords[k][:,0][::1], coords[k][:,1][::1], 10,
                       edgecolor = 'none', alpha = .25,c =ct[i],
                       label = k, )

        f.savefig(myplots.figpath( 'reinitz_exprdata_coords_nuc_offset={0}'.format(ofs)))
    if do_plot_vals:
        f = myplots.fignum(1,(8,8))
        ax = f.add_subplot(111)
        ct = mycolors.getct(len(values))
        for i,k in enumerate(values.keys()):
            ax.scatter(coords[k][:,0][::1], values[k][::1], 10,
                       edgecolor = 'none',alpha = .25,c =ct[i],
                       label = k, )

        f.savefig(myplots.figpath( 'reinitz_exprdata_ap_vals_nuc_offset={0}'.format(ofs)))

    return coords, values

import cb.p.network.bdtnp.validation as val
def get_soheil_network(max_edges = -1,
                       module_type = 'singleton'
                       ,node_restriction = None):
    mods, genes, tfs = val.get_clusters()
    tfnames, tgnames = val.fetch_cluster_mapping()
    
    mod_names = dict([(tuple(sorted([tfnames[mkey_elt -1] for mkey_elt in mkey])),
                       list(dict(m) for m in  mval))
                      for mkey, mval in mods.iteritems()])
    
    for k, v in mod_names.iteritems():
        for elt in v:
            elt.update({'gene':tgnames[elt['gene'] -1]})
            

    if module_type == 'singleton':
        mods = dict([(k[0],v) for k,v in mod_names.iteritems() if len(k) == 1])
        glam = lambda x: x['gene']
        mgenes = list(it.chain(*[[ (m,k,len(list(g))) for k,g in it.groupby(sorted(v,key = glam) ,key = glam)] 
                  for m,v in mods.iteritems() 
                  ]))
        top = sorted(mgenes, key = lambda x: x[-1])[::-1]

    
        #edge_count = dict([ (k,len(set()
        #if edge_count = -1
    else: 
        raise Exception('module type {0} not yet implemented'.format(module_type))


    if node_restriction !=  None:
        top = [e for e in top  if ( e[0] in node_restriction) and( e[1] in node_restriction)]

    if max_edges == -1: max_edges = len(top)
    edges = top[:max_edges]

    glam2 = lambda x: x[1]
    gene_groups= dict([(k,list(g)) for k,g in it.groupby(sorted(edges, key = glam2), key = glam2)])
    targets = dict([(k, {'tfs':[e[0] for e in elt], 'weights':[e[2] for e in elt]})
                    for k , elt in gene_groups.iteritems()])
    
    return nfutils.nx_from_network(targets, name = 'soheil_{0}_{1}'.format(module_type, max_edges)
                                   ,reset = True)

def check_network(net_name = 'binding', 
                  dataset_name = 'reinitz',
                  data_ofs = 4,
                  max_edges = -1):

    if dataset_name == 'reinitz':
        coords, values = get_reinitz_data(ofs = data_ofs)
    elif dataset_name == 'bdtnp':
        data = nio.getBDTNP()
        meta = nio.getBDTNP(misc = True)
        values =  dict([( k, v['vals'][:,data_ofs] ) for k,v in data.iteritems()]) 
        coords  = array([meta['x']['vals'][:,data_ofs],meta['y']['vals'][:,data_ofs]])
    
    else:
        raise Exception('data set {0} not yet implemented'.format(dataset_name))

    nets = comp.get_graphs()
    if net_name == 'binding':
        network = nets['bn']
    elif net_name =='clusters':
        network = get_soheil_network(max_edges = max_edges,
                                     node_restriction = values.keys())
    else:
        raise Exception('type not implemented: {0}'.format(net_name))

    nodes = values.keys()
    nodes_allowed = set(nodes)

    f = myplots.fignum(1,(8,8))
    ax = f.add_subplot(111)
    targets = {}

    edges = []
    
    for n in nodes:
        targets[n] = []
        if n in network:
            targets[n] = nodes_allowed.intersection(network[n].keys())
            
    xax = linspace(-1,1,20)

    edges = list(it.chain(*[[(e,v2) for v2 in v] for e, v in targets.iteritems()]))
    ccofs = [e for e in [ corrcoef(values[tf], values[tg])[0,1] for tf, tg in edges] if not isnan(e)]
    
    count, kde = make_kde(ccofs)
    

    ax.hist(ccofs,xax,label = net_name)
    h =histogram(ccofs,xax)
    ax.fill_between(xax,kde(xax)*max(h[0]),label = net_name,zorder = 1,alpha = .5)



    myplots.maketitle(ax,'edge correlations kde for {0}'.format('\n{2} data (data offset={0})\n(net_name={1})\n(max_edges={3})'
                                                                .format(data_ofs, net_name, dataset_name, max_edges) ),\
                          subtitle = 'n_edges = {0}'.format(len(edges)))
    ax.legend()
    f.savefig(myplots.figpath('network_edge_corrs_data_ofs={0}_net={1}_expr={2}_max_edges={3}'
                              .format(data_ofs,net_name,dataset_name, max_edges)))
    


def check_edges(network, nodes, data):
    pass
    
