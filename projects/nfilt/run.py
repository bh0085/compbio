import cb.p.network.io as nio
import plots as nfplots
import utils as nfu
import cb.utils.colors as mycolors

import cb.utils.plots as myplots
import matplotlib.pyplot as plt
import networkx as nx
import cb.config as cfg

import itertools as it
from numpy import *
import numpy as np

figtemplate = cfg.dataPath('figs/filter/run_{0}.pdf')


def run( reset = False,
         base_net = 'kn',
         comp_net = 'fn',
         demand_bdtnp = False):
    tgs,tfs = nio.getNet()
    ktgs,ktfs = nio.getKNet()
    bd = nio.getBDTNP()
    #btgs,btfs = nio.getBDTNP()
    sush = nio.getSush(on_fail = 'compute')
    


    tfset = set(ktfs.keys())
    tgset = set(ktgs.keys())

    
    tg_int = set(tgs.keys()).intersection(ktgs.keys())
    tf_int = set(tfs.keys()).intersection(ktfs.keys())
    
    if demand_bdtnp:
        tg_int = tg_int.intersection(bd.keys())
        tf_int = tf_int.intersection(bd.keys())



    sfRN = [(tf, tg, float(wt)) 
            for tg, elt in tgs.iteritems() if tg in tg_int
            for tf, wt  in zip(elt['tfs'], elt['weights']) if tf in tf_int]
    

    kRN = [(tf, tg, float(wt)) 
            for tg, elt in ktgs.iteritems() if tg in tg_int
            for tf, wt  in zip(elt['tfs'], elt['weights']) if tf in tf_int]

    #Sushmita network with signed edges
    suRN =  [(tf, tg, float(wt)) 
            for tg, elt in sush.iteritems() if tg in tg_int
            for tf, wt  in zip(elt['tfs'], elt['weights']) if tf in tf_int]

    #Sushmita network with unsigned edges
    suaRN = [(tf, tg, abs(float(wt))) 
             for tg, elt in sush.iteritems() if tg in tg_int
             for tf, wt  in zip(elt['tfs'], elt['weights']) if tf in tf_int]
             
    edges = [ kRN, sfRN, suRN, suaRN]


    ng = 4
    fg, kg, sug, suag = [nx.DiGraph() for i in range(4)]
    
    
    nodes = array(list(tf_int.union(tg_int)))
    graphs =  {'kg':kg,'fg':fg,'sug':sug,'suag':suag}

    for g, edges in zip(graphs.values(), edges):
        g.add_nodes_from(nodes)
        g.add_weighted_edges_from(edges)
        


    
    for gname in ['fg','suag']:
        for prc in [10,50,75,85,90,95,98,99]:
            thr = percentile([e[2]['weight'] for e in nx.to_edgelist(graphs[gname])], prc)
            graphs.update([('{0}_thr{1:2.2}'.format(gname,thr),
                            nfu.thr_graph(graphs[gname],thr))])    
            
    v0 = graphs.values()
    k0 = graphs.keys()

    tot_edges = len(nx.to_edgelist(graphs['fg']))
    for k, v in zip(k0,v0):
        for n_c in [2,4,8 ,12]:
            for max_edges in array([.5,1.,2.,5.]) * tot_edges :
                if not 'thr' in k:
                    continue
                gfilt = nfu.filter_graph(v, n_c = n_c)
                gfilt = nfu.top_edges(gfilt, max_edges = max_edges)
                gthr = nfu.thr_graph(gfilt, 1e-8)
                graphs.update([('{0}_flt{1}'.format(k,n_c),gfilt)])
                graphs.update([('{0}_flt{1}_thr0'.format(k,n_c),gthr)])
            


        

    #nfplots.show(kg,pos,node_color = 'none')
    
    #nfplots.show(fg,pos,node_color = 'white', alpha = .2, with_labels = False)
    
    return graphs


def run2( reset = False,
         base_net = 'kn',
         comp_net = 'fn',
         demand_bdtnp = False):
    bd = nio.getBDTNP()

    ktgs,ktfs = nio.getKNet()
    tgs,tfs = nio.getNet()
    sush = nio.getSush(on_fail = 'compute')
    
    tfset = set(ktfs.keys())
    tgset = set(ktgs.keys())

    
    tg_int = set(tgs.keys()).intersection(ktgs.keys())
    tf_int = set(tfs.keys()).intersection(ktfs.keys())
    
    if demand_bdtnp:
        tg_int = tg_int.intersection(bd.keys())
        tf_int = tf_int.intersection(bd.keys())
    
    if base_net =='kn':
        b_edges = [(tf, tg, float(wt)) 
               for tg, elt in ktgs.iteritems() if tg in tg_int
               for tf, wt  in zip(elt['tfs'], elt['weights']) if tf in tf_int]

    if comp_net == 'fn':
        c_edges = [(tf, tg, float(wt)) 
                for tg, elt in tgs.iteritems() if tg in tg_int
                for tf, wt  in zip(elt['tfs'], elt['weights']) if tf in tf_int]
    elif comp_net == 'sn':
        #Sushmita network with signed edges
        c_edges =  [(tf, tg, float(wt)) 
                 for tg, elt in sush.iteritems() if tg in tg_int
                 for tf, wt  in zip(elt['tfs'], elt['weights']) if tf in tf_int]
    elif comp_net == 'sna':
        #Sushmita network with unsigned edges
        c_edges = [(tf, tg, abs(float(wt))) 
                 for tg, elt in sush.iteritems() if tg in tg_int
                 for tf, wt  in zip(elt['tfs'], elt['weights']) if tf in tf_int]
    elif comp_net == 'kn':
        c_edges = [(tf, tg, float(wt)) 
           for tg, elt in ktgs.iteritems() if tg in tg_int
           for tf, wt  in zip(elt['tfs'], elt['weights']) if tf in tf_int]

    
    ng = 4    
    nodes = array(list(tf_int.union(tg_int)))

    bg = nx.DiGraph()
    bg.add_nodes_from(nodes)
    bg.add_weighted_edges_from(b_edges)

    cg = nx.DiGraph()
    cg.add_nodes_from(nodes)
    cg.add_weighted_edges_from(c_edges)

    cgraphs =  {comp_net:cg}
    v0 = cgraphs.values()
    k0 = cgraphs.keys()
   
 
    for k,g in zip(k0,v0):
        for prc in [10,50,75,85,90,95,98,99]:
            thr = percentile([e[2]['weight'] 
                              for e in nx.to_edgelist(g)], prc)
            cgraphs.update([('{0}_thr{1:2.2}'.format(k,thr),
                             nfu.thr_graph(g,thr))])    
            
    v0 = cgraphs.values()
    k0 = cgraphs.keys()

    for k, v in zip(k0,v0):
        tot_edges = len(nx.to_edgelist(v))
        for n_c in [2,4,8 ,12]:
            for max_edges in array([.5,1.,2.,5.]) * tot_edges :
                if not 'thr' in k:
                    continue
                gfilt = nfu.filter_graph(v, n_c = n_c)
                gfilt = nfu.top_edges(gfilt, max_edges = max_edges)
                gthr = nfu.thr_graph(gfilt, 1e-8)
                cgraphs.update([('{0}_flt{1}'.format(k,n_c),gfilt)])
                cgraphs.update([('{0}_flt{1}_thr0'.format(k,n_c),gthr)])
            


    '''
When you don't have hidden variables, networks can be modelled mby information criterion.

In what settings can you incur causality from datasets.

You need a prior to limit the number of arrowsin your graph:
  
The idea: come up with an idea from computational learning theory 
and come up with a model for interventions.

Spatial, genetic, time data to penalize network edges...

Granger causality uses time varying data to 
'''

        

    #nfplots.show(kg,pos,node_color = 'none')
    
    #nfplots.show(fg,pos,node_color = 'white', alpha = .2, with_labels = False)
    
    return bg, cgraphs


def gviz(graphs):
    pass
    


def run_sig(genes, show_disc = False, weighted = True):
    
    genes2 = []
    for g in genes:
        genes2.append((g[0], [[gelt[0], gelt[1]] for gelt in g[1]], g[2]))
    genes = genes2

    modules = [m[0] for m in genes]

    if len(modules[0]) == 2:
        module_type = 'doubles'
    else:
        module_type = 'triples'

    counts = [m[1] for m in genes]
    tgs,tfs = nio.getNet()
    bd = nio.getBDTNP()
    
    nodes_allowed = set(bd.keys())
    cnodes = list(nodes_allowed)

    dnodes = []
    dedges = []
    cedges = []
    cnodes = []
    for m in genes:
        for tginfo in m[1]:
            tg = tginfo[0]
            tg_mcount = tginfo[1]
            dtgnode = '{0}_{1}_mod{2}'.format(tg,tg,m[0])
            ctgnode = '{0}'.format(tg)
            dnodes.append(dtgnode)
            cnodes.append(ctgnode)
            for tf in m[0]:
                dtfnode = '{0}_{1}_mod{2}'.format(tf,tg,m[0])
                ctfnode = '{0}'.format(tf)
                dnodes.append(dtfnode)
                cnodes.append(ctfnode)
                dedges.append((dtfnode, dtgnode,tg_mcount))
                cedges.append((ctfnode, ctgnode,tg_mcount))
                
    nodes_allowed = list(set(cnodes))

    if show_disc:
        dgraph, cgraph = [nx.Graph() for i in range(2)]
        dgraph.add_nodes_from(list(set(dnodes)))
        dgraph.add_weighted_edges_from(list(set(dedges)))
        f = myplots.fignum(4, (8,8))
        ax = f.add_subplot(111)
        pos=nx.graphviz_layout(dgraph,prog="neato")
        # color nodes the same in each connected subgraph
        C=nx.connected_component_subgraphs(dgraph)
        for g in C:
            c=[random.random()]*nx.number_of_nodes(g) # random color...
            nx.draw(g,
                    pos,
                    node_size=40,
                    node_color=c,
                    vmin=0.0,
                    vmax=1.0,
                    with_labels=False
                    )
        figtitle = 'mcmc_disc'
        f.savefig(figtemplate.format(figtitle))
        return


    
    cgraph = nx.DiGraph() 
    cgraph.add_nodes_from(cnodes)
    
    cedgegrps = [(k,list(g)) for k, g in it.groupby(\
            sorted(cedges, key = lambda x: (x[0],x[1])),
            key =  lambda x: (x[0],x[1]))]
    cedges = [ (k[0],k[1], sum([gelt[2] for gelt in g])) 
                 for k,g in cedgegrps] 
                 
    if weighted == False:
        for ce in cedges:
            ce[2] = 1
        
    cgraph.add_weighted_edges_from(list(set(cedges)))


    sfRN = [(tf, tg, float(wt)) 
            for tg, elt in tgs.iteritems() if tg in nodes_allowed
            for tf, wt  in zip(elt['tfs'], elt['weights']) if tf in nodes_allowed]
    fg = nx.DiGraph()
    fg.add_nodes_from(cnodes)
    fg.add_weighted_edges_from(sfRN)


    colors = mycolors.getct(len(cnodes))

    f = myplots.fignum(5, (8,8))
    ax =f.add_subplot(111)
    pos=nx.graphviz_layout(fg,prog="neato")
    # color nodes the same in each connected subgraph
    nx.draw(cgraph,
            pos,
            node_size=100,
            node_color=colors,
            vmin=0.0,
            vmax=1.0,
            with_labels=False,
            alpha = 1.
            )
    ax.set_title('connectivity of MCMC for network {0}'.format(module_type))
    figtitle = 'mcmc_network_{0}{1}'.\
        format(module_type,'' if weighted else 'unweighted')
    f.savefig(figtemplate.format(figtitle))




    f = myplots.fignum(5, (8,8))
    ax =f.add_subplot(111)
    #pos=nx.graphviz_layout(fg,prog="neato")
    # color nodes the same in each connected subgraph
    nx.draw(fg,
            pos,
            node_size=100,
            node_color=colors,
            vmin=0.0,
            vmax=1.0,
            with_labels=False
            )


    ax.set_title('connectivity of reference for network {0}'.format(module_type))

    figtitle = 'mcmc_ref_network_{0}{1}'.\
        format(module_type,'' if weighted else 'unweighted')
    f.savefig(figtemplate.format(figtitle))

                

    graphs = {'mcmc':cgraph,'network':fg}

            
    v0 = graphs.values()
    k0 = graphs.keys()

    for k,g in zip(k0,v0):
        for prc in [1,50,95]:
            thr = percentile([e[2]['weight'] 
                              for e in nx.to_edgelist(g)], prc)
            graphs.update([('{0}_thr{1}%'.format(k,prc),
                             nfu.thr_graph(g,thr))])    
            
    v0 = graphs.values()
    k0 = graphs.keys()

    for k, v in zip(k0,v0):
        tot_edges = len(nx.to_edgelist(fg))
        for n_c in [2,4,6,8,12,20]:
            for max_edges in array([.5,1.,2.]) * tot_edges :
                gfilt = nfu.filter_graph(v, n_c = n_c)
                gfilt = nfu.top_edges(gfilt, max_edges = max_edges)
                gthr = nfu.thr_graph(gfilt, 1e-8)
                graphs.update([('{0}_flt{1}'.format(k,n_c),gfilt)])
                graphs.update([('{0}_flt{1}_thr0'.format(k,n_c),gthr)])
            
    return graphs

def show_mcmc(graphs,nc_base =12,
              weighted = True, module_type = 'doubles',
              measure = 'cosine'):

    nodelist = graphs.values()[0].nodes()

    print 'using module type: {0}'.format(module_type)
    adjs = [ array(nx.adj_matrix(g, nodelist = nodelist)) for g in graphs.values() ]
    nrms = []
    for a in adjs:
            n = sqrt(sum(a**2))
            nrms.append(a / n)
    

    idxs_allowed =[ i for i, k in  enumerate(graphs.keys()) if 'mcmc' in k]
    
    belt = graphs.keys().index('network_flt{0}'.format(nc_base))
    if measure == 'cosine':
        sims = array([round(nfu.cosine_adj(a1,nrms[belt]),8) 
                      for i, a1 in enumerate(nrms) if i in idxs_allowed])
    elif measure == 'dotprod':
        sims = array([nfu.dotprod(a1,adjs[belt]) 
                      for i, a1 in enumerate(adfs) if i in idxs_allowed])
    else:
        raise Exception()

    keys_allowed = [k for i, k in enumerate(graphs.keys()) 
                    if i in idxs_allowed]
    srto = argsort([k for i, k in enumerate(graphs.keys()) 
                    if i in idxs_allowed]) 
    #XVALs give ranks of each key index.
    xvals = argsort(srto)


    cols = map(lambda x: 
               ('flt' in x and x.count('thr') > 1) and 'orange' or
               ('flt' in x) and 'red' or
               ('thr' in x) and 'yellow' or
               ('fg' in x) and 'green' or 
               ('su' in x) and 'blue' or 
               'black', keys_allowed)

    yvals = sims

    f = plt.gcf()
    f = myplots.fignum(3, (.25 * len(sims),10))
    f.clear()
    ax = f.add_subplot(111)
    myplots.padded_limits(ax,xvals,yvals + [0.], margin = [.02,.02])
    ax.scatter(xvals,yvals,100, color = cols)
    ax.set_ylabel('red fly similarity ({0})'.format(measure))
    ax.set_xlabel('networks')
    ax.set_xticklabels([])
    ax.set_xticks([])

    mark_ys = [0, median(sims), mean(sims), sort(sims)[::-1][1],1]
    ax.hlines(mark_ys, *ax.get_xlim(), linestyle = ':',alpha = .2)
    ax.vlines(range(len(xvals))[::10], *ax.get_ylim(), linestyle = ':',alpha = .1)
    

    

    figtitle = 'mcmc_net_comparisons_{0}_{1}net_comps{2}_sim_{3}_nolabels'.\
        format(module_type,nc_base ,'' if weighted else 'unweighted',measure)
    f.savefig(figtemplate.format(figtitle))

                

    ax.set_xticks(range(len(srto)))
    #ax.annotate('\n'.join(' '.join(z) for z in zip(graphs.keys())),
    #            [0,1],xycoords = 'axes fraction', va = 'top')
    
    ax.set_xticklabels([keys_allowed[i] for i in srto], 
                       rotation = 90, size = 'xx-small',
                       )#color = cols)


    figtitle = 'mcmc_net_comparisons_{0}_{1}net_comps{2}_sim_{3}_labels'.\
        format(module_type,nc_base ,'' if weighted else 'unweighted', measure)
    f.savefig(figtemplate.format(figtitle))


    
import cb.p.network.io as nio
import plots as nfplots
import utils as nfu
import cb.utils.colors as mycolors

import cb.utils.plots as myplots
import matplotlib.pyplot as plt
import networkx as nx
import cb.config as cfg

import itertools as it
from numpy import *
import numpy as np

figtemplate = cfg.dataPath('figs/filter/run_{0}.pdf')


def run( reset = False,
         base_net = 'kn',
         comp_net = 'fn',
         demand_bdtnp = False):
    tgs,tfs = nio.getNet()
    ktgs,ktfs = nio.getKNet()
    bd = nio.getBDTNP()
    #btgs,btfs = nio.getBDTNP()
    sush = nio.getSush(on_fail = 'compute')
    


    tfset = set(ktfs.keys())
    tgset = set(ktgs.keys())

    
    tg_int = set(tgs.keys()).intersection(ktgs.keys())
    tf_int = set(tfs.keys()).intersection(ktfs.keys())
    
    if demand_bdtnp:
        tg_int = tg_int.intersection(bd.keys())
        tf_int = tf_int.intersection(bd.keys())



    sfRN = [(tf, tg, float(wt)) 
            for tg, elt in tgs.iteritems() if tg in tg_int
            for tf, wt  in zip(elt['tfs'], elt['weights']) if tf in tf_int]
    

    kRN = [(tf, tg, float(wt)) 
            for tg, elt in ktgs.iteritems() if tg in tg_int
            for tf, wt  in zip(elt['tfs'], elt['weights']) if tf in tf_int]

    #Sushmita network with signed edges
    suRN =  [(tf, tg, float(wt)) 
            for tg, elt in sush.iteritems() if tg in tg_int
            for tf, wt  in zip(elt['tfs'], elt['weights']) if tf in tf_int]

    #Sushmita network with unsigned edges
    suaRN = [(tf, tg, abs(float(wt))) 
             for tg, elt in sush.iteritems() if tg in tg_int
             for tf, wt  in zip(elt['tfs'], elt['weights']) if tf in tf_int]
             
    edges = [ kRN, sfRN, suRN, suaRN]


    ng = 4
    fg, kg, sug, suag = [nx.DiGraph() for i in range(4)]
    
    
    nodes = array(list(tf_int.union(tg_int)))
    graphs =  {'kg':kg,'fg':fg,'sug':sug,'suag':suag}

    for g, edges in zip(graphs.values(), edges):
        g.add_nodes_from(nodes)
        g.add_weighted_edges_from(edges)
        


    
    for gname in ['fg','suag']:
        for prc in [10,50,75,85,90,95,98,99]:
            thr = percentile([e[2]['weight'] for e in nx.to_edgelist(graphs[gname])], prc)
            graphs.update([('{0}_thr{1:2.2}'.format(gname,thr),
                            nfu.thr_graph(graphs[gname],thr))])    
            
    v0 = graphs.values()
    k0 = graphs.keys()

    tot_edges = len(nx.to_edgelist(graphs['fg']))
    for k, v in zip(k0,v0):
        for n_c in [2,4,8 ,12]:
            for max_edges in array([.5,1.,2.,5.]) * tot_edges :
                if not 'thr' in k:
                    continue
                gfilt = nfu.filter_graph(v, n_c = n_c)
                gfilt = nfu.top_edges(gfilt, max_edges = max_edges)
                gthr = nfu.thr_graph(gfilt, 1e-8)
                graphs.update([('{0}_flt{1}'.format(k,n_c),gfilt)])
                graphs.update([('{0}_flt{1}_thr0'.format(k,n_c),gthr)])
            


        

    #nfplots.show(kg,pos,node_color = 'none')
    
    #nfplots.show(fg,pos,node_color = 'white', alpha = .2, with_labels = False)
    
    return graphs


def run2( reset = False,
         base_net = 'kn',
         comp_net = 'fn',
         demand_bdtnp = False):
    bd = nio.getBDTNP()

    ktgs,ktfs = nio.getKNet()
    tgs,tfs = nio.getNet()
    sush = nio.getSush(on_fail = 'compute')
    
    tfset = set(ktfs.keys())
    tgset = set(ktgs.keys())

    
    tg_int = set(tgs.keys()).intersection(ktgs.keys())
    tf_int = set(tfs.keys()).intersection(ktfs.keys())
    
    if demand_bdtnp:
        tg_int = tg_int.intersection(bd.keys())
        tf_int = tf_int.intersection(bd.keys())
    
    if base_net =='kn':
        b_edges = [(tf, tg, float(wt)) 
               for tg, elt in ktgs.iteritems() if tg in tg_int
               for tf, wt  in zip(elt['tfs'], elt['weights']) if tf in tf_int]

    if comp_net == 'fn':
        c_edges = [(tf, tg, float(wt)) 
                for tg, elt in tgs.iteritems() if tg in tg_int
                for tf, wt  in zip(elt['tfs'], elt['weights']) if tf in tf_int]
    elif comp_net == 'sn':
        #Sushmita network with signed edges
        c_edges =  [(tf, tg, float(wt)) 
                 for tg, elt in sush.iteritems() if tg in tg_int
                 for tf, wt  in zip(elt['tfs'], elt['weights']) if tf in tf_int]
    elif comp_net == 'sna':
        #Sushmita network with unsigned edges
        c_edges = [(tf, tg, abs(float(wt))) 
                 for tg, elt in sush.iteritems() if tg in tg_int
                 for tf, wt  in zip(elt['tfs'], elt['weights']) if tf in tf_int]
    elif comp_net == 'kn':
        c_edges = [(tf, tg, float(wt)) 
           for tg, elt in ktgs.iteritems() if tg in tg_int
           for tf, wt  in zip(elt['tfs'], elt['weights']) if tf in tf_int]

    
    ng = 4    
    nodes = array(list(tf_int.union(tg_int)))

    bg = nx.DiGraph()
    bg.add_nodes_from(nodes)
    bg.add_weighted_edges_from(b_edges)

    cg = nx.DiGraph()
    cg.add_nodes_from(nodes)
    cg.add_weighted_edges_from(c_edges)

    cgraphs =  {comp_net:cg}
    v0 = cgraphs.values()
    k0 = cgraphs.keys()
   
 
    for k,g in zip(k0,v0):
        for prc in [10]:
            thr = percentile([e[2]['weight'] 
                              for e in nx.to_edgelist(g)], prc)
            cgraphs.update([('{0}_thr{1:2.2}'.format(k,thr),
                             nfu.thr_graph(g,thr))])    
            gt =  nfu.thr_graph(g,thr)
    v0 = cgraphs.values()
    k0 = cgraphs.keys()

    for k, v in zip(k0,v0):
        tot_edges = len(nx.to_edgelist(v))
        for n_c in [1,2,4]:
            for max_edges in array([.2,.5,1.]) * tot_edges :
                if not 'thr' in k:
                    continue
                gfilt = nfu.filter_graph(v, n_c = n_c)
                gfilt = nfu.top_edges(gfilt, max_edges = max_edges)
                gthr = nfu.thr_graph(gfilt, 1e-8)
                cgraphs.update([('{0}_flt{1}_nedge{2}'.format(k,n_c,max_edges),gfilt)])
                cgraphs.update([('{0}_flt{1}_nedge{2}_thr0'.format(k,n_c,max_edges),gthr)])

    '''
When you don't have hidden variables, networks can be modelled mby information criterion.

In what settings can you incur causality from datasets.

You need a prior to limit the number of arrowsin your graph:
  
The idea: come up with an idea from computational learning theory 
and come up with a model for interventions.

Spatial, genetic, time data to penalize network edges...

Granger causality uses time varying data to 
'''

        

    #nfplots.show(kg,pos,node_color = 'none')
    
    #nfplots.show(fg,pos,node_color = 'white', alpha = .2, with_labels = False)
    
    return bg, cgraphs


def gviz(graphs):
    pass
    
def gdraw(bgraph,cgraphs, plotname = 'default_name', measure = 'cosine'):
    #pos = nx.graphviz_layout(cgraphs['kg'])

    nodelist = bgraph.nodes()

    adjs = [ array(nx.adj_matrix(g, nodelist = nodelist)) for g in cgraphs.values() ]



    badj =  array(nx.adj_matrix(bgraph, nodelist = nodelist))
    bnrm = badj / sqrt(sum(badj**2))

    if measure == 'cosine':
        nrms = []
        bnrm = badj / sqrt(sum(badj**2))
        for a in adjs:
            n = sqrt(sum(a**2))
            nrms.append(a / n)
    
        sims = array([round(nfu.cosine_adj(a1,bnrm),8) for a1 in nrms])
    elif measure =='jaccard':
        sims = array([round(nfu.dotprod(a1,badj),8)/ (sum(a1) + sum(badj)) 
                      for a1 in adjs])

    elif measure =='specificity':
        sims = array([round(nfu.dotprod(a1,badj),8)/sum(a1) for a1 in adjs])

    elif measure =='sensitivity':
        sims = array([round(nfu.dotprod(a1,badj),8)/sum(badj) for a1 in adjs])

    else:
        raise Exception()




    srto = argsort(cgraphs.keys()) 
    #XVALs give ranks of each key index.
    xvals = argsort(srto)


    cols = map(lambda x: 
               ('flt' in x and x.count('thr') > 1) and 'orange' or
               ('flt' in x) and 'red' or
               ('thr' in x) and 'yellow' or
               ('fg' in x) and 'green' or 
               ('su' in x) and 'blue' or 
               'black', cgraphs.keys())

    yvals = sims

    f = plt.gcf()
    f = myplots.fignum(3, (.25 * len(sims),10))
    f.clear()
    ax = f.add_subplot(111)
    myplots.padded_limits(ax,xvals,yvals + [0.], margin = [.02,.02])
    ax.scatter(xvals,yvals,100, color = cols)
    ax.set_ylabel('red fly similarity ({0})'.format(measure))
    ax.set_xlabel('networks')
    ax.set_xticklabels([])
    ax.set_xticks([])

    mark_ys = [0, median(sims), mean(sims), sort(sims)[::-1][1],1]
    ax.hlines(mark_ys, *ax.get_xlim(), linestyle = ':',alpha = .2)
    

    f.savefig(cfg.dataPath('figs/meAWG/filter_{0}_meth_{1}_nolabels.pdf'.\
                               format(plotname,measure)))


    ax.set_xticks(range(len(srto)))
    cols_added = []
    annotes = []
    for z in zip(cgraphs.keys(),cols):
        if not z[1] in cols_added:
            annotes.append( ' '.join(z))
            cols_added.append(z[1])
        
    ax.annotate('\n'.join(annotes), 
                [1,1],xycoords = 'axes fraction', va = 'top', ha = 'right')
    
    ax.set_xticklabels([cgraphs.keys()[i] for i in srto], 
                       rotation = 90, va = 'bottom', size = 'xx-small')

    f.savefig(cfg.dataPath('figs/meAWG/filter_{0}_meth_{1}_labels.pdf'.\
                               format(plotname,measure)))

def gdraw0(graphs, plotname = 'default_name', measure = 'cosine'):
    pos = nx.graphviz_layout(graphs['kg'])


    adjs = [ array(nx.adj_matrix(g)) for g in graphs.values() ]
    nrms = []
    for a in adjs:
            n = sqrt(sum(a**2))
            nrms.append(a / n)
    
    kgelt = graphs.keys().index('kg')
    if measure == 'cosine':
        sims = array([round(nfu.cosine_adj(a1,nrms[kgelt]),8) for a1 in nrms])
    else:
        raise Exception()

    kg = graphs['kg']
    srto = argsort(graphs.keys()) 
    #XVALs give ranks of each key index.
    xvals = argsort(srto)


    cols = map(lambda x: 
               ('flt' in x and x.count('thr') > 1) and 'orange' or
               ('flt' in x) and 'red' or
               ('thr' in x) and 'yellow' or
               ('fg' in x) and 'green' or 
               ('su' in x) and 'blue' or 
               'black', graphs.keys())

    yvals = sims

    f = plt.gcf()
    f = myplots.fignum(3, (.25 * len(sims),10))
    f.clear()
    ax = f.add_subplot(111)
    myplots.padded_limits(ax,xvals,yvals + [0.], margin = [.02,.02])
    ax.scatter(xvals,yvals,100, color = cols)
    ax.set_ylabel('red fly similarity ({0})'.format(measure))
    ax.set_xlabel('networks')
    ax.set_xticklabels([])
    ax.set_xticks([])
    mark_ys = [0, median(sims), mean(sims), sort(sims)[::-1][1],1]
    ax.hlines(mark_ys, *ax.get_xlim(), linestyle = ':',alpha = .2)
    

    f.savefig(cfg.dataPath('figs/meAWG/filter_{0}_meth_{1}_nolabels.pdf'.\
                               format(plotname,measure)))


    ax.set_xticks(range(len(srto)))
    ax.annotate('\n'.join(' '.join(z) for z in zip(graphs.keys(),cols)),
                [0,1],xycoords = 'axes fraction', va = 'top')
    
    ax.set_xticklabels([graphs.keys()[i] for i in srto], 
                       rotation = 45, size = 'xx-small',ha = 'right')

    f.savefig(cfg.dataPath('figs/meAWG/filter_{0}_meth_{1}_labels.pdf'.\
                               format(plotname,measure)))


def run_sig(genes, show_disc = False, weighted = True):
    
    genes2 = []
    for g in genes:
        genes2.append((g[0], [[gelt[0], gelt[1]] for gelt in g[1]], g[2]))
    genes = genes2

    modules = [m[0] for m in genes]

    if len(modules[0]) == 2:
        module_type = 'doubles'
    else:
        module_type = 'triples'

    counts = [m[1] for m in genes]
    tgs,tfs = nio.getNet()
    bd = nio.getBDTNP()
    
    nodes_allowed = set(bd.keys())
    cnodes = list(nodes_allowed)

    dnodes = []
    dedges = []
    cedges = []
    cnodes = []
    for m in genes:
        for tginfo in m[1]:
            tg = tginfo[0]
            tg_mcount = tginfo[1]
            dtgnode = '{0}_{1}_mod{2}'.format(tg,tg,m[0])
            ctgnode = '{0}'.format(tg)
            dnodes.append(dtgnode)
            cnodes.append(ctgnode)
            for tf in m[0]:
                dtfnode = '{0}_{1}_mod{2}'.format(tf,tg,m[0])
                ctfnode = '{0}'.format(tf)
                dnodes.append(dtfnode)
                cnodes.append(ctfnode)
                dedges.append((dtfnode, dtgnode,tg_mcount))
                cedges.append((ctfnode, ctgnode,tg_mcount))
                
    nodes_allowed = list(set(cnodes))

    if show_disc:
        dgraph, cgraph = [nx.Graph() for i in range(2)]
        dgraph.add_nodes_from(list(set(dnodes)))
        dgraph.add_weighted_edges_from(list(set(dedges)))
        f = myplots.fignum(4, (8,8))
        ax = f.add_subplot(111)
        pos=nx.graphviz_layout(dgraph,prog="neato")
        # color nodes the same in each connected subgraph
        C=nx.connected_component_subgraphs(dgraph)
        for g in C:
            c=[random.random()]*nx.number_of_nodes(g) # random color...
            nx.draw(g,
                    pos,
                    node_size=40,
                    node_color=c,
                    vmin=0.0,
                    vmax=1.0,
                    with_labels=False
                    )
        figtitle = 'mcmc_disc'
        f.savefig(figtemplate.format(figtitle))
        return


    
    cgraph = nx.DiGraph() 
    cgraph.add_nodes_from(cnodes)
    
    cedgegrps = [(k,list(g)) for k, g in it.groupby(\
            sorted(cedges, key = lambda x: (x[0],x[1])),
            key =  lambda x: (x[0],x[1]))]
    cedges = [ (k[0],k[1], sum([gelt[2] for gelt in g])) 
                 for k,g in cedgegrps] 
                 
    if weighted == False:
        for ce in cedges:
            ce[2] = 1
        
    cgraph.add_weighted_edges_from(list(set(cedges)))


    sfRN = [(tf, tg, float(wt)) 
            for tg, elt in tgs.iteritems() if tg in nodes_allowed
            for tf, wt  in zip(elt['tfs'], elt['weights']) if tf in nodes_allowed]
    fg = nx.DiGraph()
    fg.add_nodes_from(cnodes)
    fg.add_weighted_edges_from(sfRN)


    colors = mycolors.getct(len(cnodes))

    f = myplots.fignum(5, (8,8))
    ax =f.add_subplot(111)
    pos=nx.graphviz_layout(fg,prog="neato")
    # color nodes the same in each connected subgraph
    nx.draw(cgraph,
            pos,
            node_size=100,
            node_color=colors,
            vmin=0.0,
            vmax=1.0,
            with_labels=False,
            alpha = 1.
            )
    ax.set_title('connectivity of MCMC for network {0}'.format(module_type))
    figtitle = 'mcmc_network_{0}{1}'.\
        format(module_type,'' if weighted else 'unweighted')
    f.savefig(figtemplate.format(figtitle))




    f = myplots.fignum(5, (8,8))
    ax =f.add_subplot(111)
    #pos=nx.graphviz_layout(fg,prog="neato")
    # color nodes the same in each connected subgraph
    nx.draw(fg,
            pos,
            node_size=100,
            node_color=colors,
            vmin=0.0,
            vmax=1.0,
            with_labels=False
            )


    ax.set_title('connectivity of reference for network {0}'.format(module_type))

    figtitle = 'mcmc_ref_network_{0}{1}'.\
        format(module_type,'' if weighted else 'unweighted')
    f.savefig(figtemplate.format(figtitle))

                

    graphs = {'mcmc':cgraph,'network':fg}

            
    v0 = graphs.values()
    k0 = graphs.keys()

    for k,g in zip(k0,v0):
        for prc in [1,50,95]:
            thr = percentile([e[2]['weight'] 
                              for e in nx.to_edgelist(g)], prc)
            graphs.update([('{0}_thr{1}%'.format(k,prc),
                             nfu.thr_graph(g,thr))])    
            
    v0 = graphs.values()
    k0 = graphs.keys()

    for k, v in zip(k0,v0):
        tot_edges = len(nx.to_edgelist(fg))
        for n_c in [2,4,6,8,12,20]:
            for max_edges in array([.5,1.,2.]) * tot_edges :
                gfilt = nfu.filter_graph(v, n_c = n_c)
                gfilt = nfu.top_edges(gfilt, max_edges = max_edges)
                gthr = nfu.thr_graph(gfilt, 1e-8)
                graphs.update([('{0}_flt{1}'.format(k,n_c),gfilt)])
                graphs.update([('{0}_flt{1}_thr0'.format(k,n_c),gthr)])
            
    return graphs

def show_mcmc(graphs,nc_base =12,
              weighted = True, module_type = 'doubles',
              measure = 'cosine'):

    print 'using module type: {0}'.format(module_type)
    adjs = [ array(nx.adj_matrix(g)) for g in graphs.values() ]
    nrms = []
    for a in adjs:
            n = sqrt(sum(a**2))
            nrms.append(a / n)
    

    idxs_allowed =[ i for i, k in  enumerate(graphs.keys()) if 'mcmc' in k]
    
    belt = graphs.keys().index('network_flt{0}'.format(nc_base))
    if measure == 'cosine':
        sims = array([round(nfu.cosine_adj(a1,nrms[belt]),8) 
                      for i, a1 in enumerate(nrms) if i in idxs_allowed])
    elif measure == 'dotprod':
        sims = array([nfu.dotprod(a1,adjs[belt]) 
                      for i, a1 in enumerate(adfs) if i in idxs_allowed])
    else:
        raise Exception()

    keys_allowed = [k for i, k in enumerate(graphs.keys()) 
                    if i in idxs_allowed]
    srto = argsort([k for i, k in enumerate(graphs.keys()) 
                    if i in idxs_allowed]) 
    #XVALs give ranks of each key index.
    xvals = argsort(srto)


    cols = map(lambda x: 
               ('flt' in x and x.count('thr') > 1) and 'orange' or
               ('flt' in x) and 'red' or
               ('thr' in x) and 'yellow' or
               ('fg' in x) and 'green' or 
               ('su' in x) and 'blue' or 
               'black', keys_allowed)

    yvals = sims

    f = plt.gcf()
    f = myplots.fignum(3, (.25 * len(sims),10))
    f.clear()
    ax = f.add_subplot(111)
    myplots.padded_limits(ax,xvals,yvals + [0.], margin = [.02,.02])
    ax.scatter(xvals,yvals,100, color = cols)
    ax.set_ylabel('red fly similarity ({0})'.format(measure))
    ax.set_xlabel('networks')
    ax.set_xticklabels([])
    ax.set_xticks([])

    mark_ys = [0, median(sims), mean(sims), sort(sims)[::-1][1],1]
    ax.hlines(mark_ys, *ax.get_xlim(), linestyle = ':',alpha = .2)
    ax.vlines(range(len(xvals))[::10], *ax.get_ylim(), linestyle = ':',alpha = .1)
    

    

    figtitle = 'mcmc_net_comparisons_{0}_{1}net_comps{2}_sim_{3}_nolabels'.\
        format(module_type,nc_base ,'' if weighted else 'unweighted',measure)
    f.savefig(figtemplate.format(figtitle))

                

    ax.set_xticks(range(len(srto)))
    #ax.annotate('\n'.join(' '.join(z) for z in zip(graphs.keys())),
    #            [0,1],xycoords = 'axes fraction', va = 'top')
    
    ax.set_xticklabels([keys_allowed[i] for i in srto], 
                       rotation = 45, size = 'xx-small',
                       )#color = cols)


    figtitle = 'mcmc_net_comparisons_{0}_{1}net_comps{2}_sim_{3}_labels'.\
        format(module_type,nc_base ,'' if weighted else 'unweighted', measure)
    f.savefig(figtemplate.format(figtitle))


    
