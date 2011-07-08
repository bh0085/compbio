import networkx as nx
from numpy import *
import numpy as np
import cb.p.network.io as nio
import cb.utils.plots as myplots
import cb.utils.gdraw as gd
import cb.config as cfg

import itertools as it
figtemplate = cfg.dataPath('figs/soheil/module_signs_{0}.pdf')

def view0(modules,
          data_src = 'bdtnp',
          net_src = 'fRN',
          max_rank = 4,
          module_type = 'doubles'):
    '''
    A routine to view the sign of interaction coefficients for 
    a given transcription factor split per-cluster and per-module
    size.

    Designed to be run on the output of view_output.modules()

'''
    #COMPUTE BULK STATISTICS FOR EACH TF
    bd_data = nio.getBDTNP()
    genes = bd_data.keys()

    tfs =sorted(set(it.chain(*[k for k in modules.keys()])))
    tf_net = nx.Graph()
    tf_net.add_nodes_from(tfs)
    tf_edges = it.chain(*[[(e0, e1)
                for e0 in term
                for e1 in term
                if e0 != e1]
                for term in modules.keys()])
    tf_net.add_edges_from(tf_edges)
 
    pos = nx.graphviz_layout(tf_net)
    fig = myplots.fignum(1, (8,8))


    tfnodes = tf_net.nodes()
    tfnames = tfnodes

    all_coefs = \
        list(it.chain(*[t['coefs'] for t in modules.values()]))
    cstd = std(all_coefs)
    def colorfun(coefs):
        coef = median(coefs)
        arr = array([1.,0.,0.]) if coef < 0\
            else array([0.,0.,1.])
        return arr*min([1, abs(coef/cstd*2)])

    def widthfun(coefs, maxwid):
        return max([1,5.* len(coefs) /maxwid])

        
  

    for tf in tfs:
        fig.clf()
        ax = fig.add_subplot(111)

        ecols,ewids, estyles, ealphas =\
            [{} for i in range(4)]
        edges = []
        tf_doublet_terms = [(k, v)
                            for k, v in modules.iteritems()
                            if tf in k and len(set(k)) == 2 ] 

        tf_triplet_terms = [(k, v)
                            for k, v in modules.iteritems()
                            if tf in k and len(set(k)) == 3 ] 
        
        isdub = dict([(k,0.) for k in tfnames])
        istrip =dict([(k,0.) for k in tfnames])

        coflens = [len(e[1]['coefs']) 
                   for e in  tf_triplet_terms + tf_doublet_terms]
        max_l = max(coflens)

        for tt,v in sorted(tf_triplet_terms,
                         key = lambda x: len(x[1]['genes']))[::-1][:max_rank]:
            partners =tuple( [t for t in set(tt) if not t==tf])
            for p in partners: istrip[p] = 1 #istrip[[tfnodes.index(p) for p in partners]] = 1
            edge = partners 
            ecols[edge] = colorfun(modules[tt]['coefs'])
            ewids[edge] = widthfun(modules[tt]['coefs'],max_l)
            ealphas[edge] = 1 if module_type in ['triples','all'] \
                else .1
            estyles[edge] = 'dotted'
            
            edges.append(edge)
            
                        
        for td,v in sorted(tf_doublet_terms,
                         key = lambda x: len(x[1]['genes']))[::-1][:max_rank]:
            partners = tuple([t for t in set(td) if not t==tf])
            for p in partners: isdub[p] = 1
            edge = (tuple([tf] + list(partners))) 
            ecols[edge] = colorfun(modules[td]['coefs'])
            ewids[edge] = widthfun(modules[td]['coefs'], max_l)
            ealphas[edge] = 1 if module_type in ['doubles','all'] \
                else .1

            estyles[edge] = 'solid'

            edges.append(edge)
        
        
        tf_graph = nx.DiGraph()
        tf_graph.add_nodes_from(tfnodes)
        tf_graph.add_edges_from(edges)
        
        ckw = dict([ (k, dict(color = ecols[k], #array([1,0,0])*isdub[k] +\
                                 # array([0,0,1])*istrip[k],
                              alpha = ealphas[k],
                              linestyle = estyles[k],
                              linewidth = ewids[k],
                              arrowstyle = '-'))
                       for k in tf_graph.edges()])        
        circlepoints = dict([ (k, dict(facecolor ='white', #array([1,0,0])*isdub[k] +\
                                       # array([0,0,1])*istrip[k],
                                       alpha = round(ealphas[k],2),
                                       edgecolor = ecols[k],
                                       linestyle = estyles[k],
                                       linewidth = 3,))
                       for k in tf_graph.edges()])


        ax.set_title('Top modules for TF: {0}'.format(tf))
        myplots.padded_limits(ax,  *zip(*pos.values()))
        nodes = tf_graph.nodes()
        gd.draw(tf_graph,pos,edges, 
                scatter_nodes =tf_graph.nodes(),
                skw = {'s':[200 if n == tf else 2
                            for n in nodes],
                       'facecolor':[colorfun(modules[tuple([n])]['coefs']) 
                                if n == tf
                                else 'black' 
                                for n in nodes],
                       'edgecolor':'black',
                       'linewidth':2,
                       'alpha':1},#[1 if n == tf else .1 for n in nodes]},
                ckw = {})
        colors,alphas,tripoints = [{} for i in range(3)]
        for e in edges:
            colors[e] = 'black'
            alphas[e] = .2
            tripoints[e] = array(pos[tf])
        
        plot_tri = False
        plot_circ = True
        if plot_tri:
            gd.overlay(tf_graph, pos, edges,
                       tripoints = tripoints,
                       colors = colors,
                       alphas = alphas)
        if plot_circ:
            gd.overlay(tf_graph, pos, edges,
                       circlepoints =circlepoints)

        ax2 = fig.add_axes([.05,.05,.2,.2])
        ax2.set_xticks([])
        ax2.set_yticks([])
        coefs = modules[tuple([tf])]['coefs']
        l = len(coefs)
        sx = sy = ceil(sqrt(l))
        xs = mod(arange(l), sx)
        ys = floor( arange(l) / sx)
        cs = [colorfun([c/2]) for c in sorted(coefs)]
        ss = 100

        ax2.scatter(xs, ys, s = ss, color = cs)

        fig.savefig(figtemplate.format('tf_{0}_net_rank{1}_{2}'.\
                                           format(tf,max_rank,module_type)))

