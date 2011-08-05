'''
A few utilities to compare different networks to eachother.

Using this package, I want to be able to:
(1) compute and plot pure netowrk statistics such as clustering cofs,
transitivity,...

(2) compute and plot biological statistics such as expression correlations
averaged over network edges, GO term clustering, spatial clustering on the
BDTNP embryos.

(3) compute the same statistics for filtered versions of input networks to
ask: how do networks change in both their graph theoretic characteristics
and their biological characteristics after my filters have been applied.

'''


import cb.projects.network.io as nio
import numpy as np
from numpy import *
import cb.utils.plots as myplots
import cb.config as cfg
import cb.utils.colors as mycolors
import cb.p.nfilt.utils as nfu
import networkx as nx
import matplotlib.pyplot as plt

from scipy.stats import linregress

ecount = 200000
run_name = 'bn_edges'
figtemplate = cfg.dataPath('figs/filtering/comparator-{0}-{{0}}.pdf'.\
                               format(run_name))
net_types = 'all'
pctype = 'bn'



def get_graphs(selector = 'just_nets',
               net_choice = net_types,
               pctype = pctype,
               fixed_ecount = ecount):

    '''Get a list of graphs.

net_choices:
  'all': fetch all available nets and build graphs according to 
         selectors
  'unsup': just fetch the unsup net.

selectors:
  'just_nets': return untransformed graphs of each kind.
  'pcs_fix_ecount': return filtered graphs with the same
                   number of edges as the original input.
                   [kwds]: pctype - which graph to compute the pc of

'''



    if net_choice == 'all':
        nets  = {'unsup':nio.getNet(net_name = 'unsup')[0],
                 'logistic':nio.getNet(net_name = 'logistic')[0],
                 'kn':nio.getKNet()[0],
                 'bn':nio.getBNet()[0],
                 'mn':nio.getMNet()[0]}
    elif net_choice == 'unsup':
        nets = {'unsup':nio.getNet(net_name = 'unsup')[0]}
    elif net_choice in ['kn','bn','mn']:
        if net_choice == 'kn': nets = {'kn':nio.getKNet()[0]}
        if net_choice == 'mn': nets = {'mn':nio.getMNet()[0]} 
        if net_choice == 'bn': nets = {'bn':nio.getBNet()[0]} 
    else:
        raise Exception()

    if selector == 'just_nets':
        graphs = dict([(k, nfu.nx_from_network(v, name = k))
                       for k, v in nets.iteritems()])

    elif selector == 'pcs_fix_ecount':
        '''Compute the PC projection with an edge count fixed at the size of
        the input network.'''
        
        graphs = dict([(k, nfu.top_edges(\
                        nfu.nx_from_network(v, name = k, reset = True),
                        fixed_ecount))
                       for k, v in nets.iteritems()])
        
        for n_c in range(1, 13,3):
            graphs['{0}_nc={1}'.format(pctype,n_c)]=\
                nfu.filter_sparse(graphs[pctype],n_c=n_c, 
                                  max_edges = len(graphs[pctype].edges()),
                                  last_component = False)
            nets['{0}_nc={1}'.format(pctype,n_c)] = nets[pctype]
    else: 
        raise Exception()
    

    return graphs 
    


def compare(graphs,
            do_transitivity = False,
            do_corrs = False, 
            do_clustering = True,
            do_edge_comparisons = True,
            expr_type = 'time_course',
            edge_restriction = 'none',
            graph_selector = 'just_nets',
            bgraphs = None):
    '''
    Compare the different nets from binding, motifs, knowledge and modencode.
    '''


    if do_transitivity: make_transitivity(graphs)
    if do_corrs: make_corrplot(graphs, 
                               expr_type = expr_type,
                               edge_restriction = edge_restriction,
                               selector = graph_selector
                               )
    if do_clustering: make_clustering(graphs)
    if do_edge_comparisons: make_edge_comparisons(graphs, bgraphs)

def make_corrplot(graphs, expr_type = 'bdtnp', 
                  verbose = False,
                  edge_restriction = 'none'):
    f = myplots.fignum(3,(8,8))
    
    if expr_type == 'time_course':
        xpr = nio.getTC()
    elif expr_type=='cell_line':
        xpr = nio.getCL()
    elif expr_type=='bdtnp':
        xpr = nio.getBDTNP()
    else:
        raise Exception()
    
    if edge_restriction == 'bdtnp':
        edges_allowed = set(nio.getBDTNP().keys())

    kdes = {}
    counts = {}

    for k, v in graphs.iteritems():
        vals = []
        for tg, tf in v.edges()        :
            if edge_restriction == 'none': 
                pass
            elif edge_restriction == 'bdtnp':
                if not tg in edges_allowed or\
                        not tf in edges_allowed:
                    continue

            if not xpr.has_key(tf) or not xpr.has_key(tg): continue
            tfx = xpr[tf]
            tgx = xpr[tg]
            #Deal with the BDTNP data that may lack certain points:
            if type(tfx) == dict:
                times_allowed = \
                    nonzero(np.min(array([sum(tfx['vals'], 0),\
                                              sum(tgx['vals'],0)]),0))[0]
                tgx = reshape(tgx['vals'][times_allowed, :],-1)
                tfx = reshape(tfx['vals'][times_allowed, :],-1)

            if len(tgx) == 0: continue
            if not max(tgx) - min(tgx) or not max(tfx) - min(tfx):
                if verbose: print 'TGX or TFX unchanging... SKIPPING {0}'.format((tfx,tgx))
                continue

            #lin = linregress(tfx, tgx)
            #rsquared = lin[2]**2
            ccof = corrcoef(tfx, tgx)[0,1]
        
            
            if not isfinite(ccof): 
                print 'problem for {0}'.format((tf,tg))
                continue
            #print 'sxs'

            vals.append(ccof)

                    
        from scipy.stats import gaussian_kde
        counts[k] = len(vals)
        dist = array(vals)
        density = gaussian_kde(dist)
        density.covariance_factor = lambda : .25
        density._compute_covariance()
        kdes[k] = (density)


    ax = f.add_subplot(111)
    plots = []
    labels = []
    for k, kde  in kdes.iteritems():
        xax = linspace(-1,1,100)
        plots.append(ax.plot(xax, kde(xax), 
                             label  = k)[0])
        labels.append('{0}, $n_e={1}$'.format(k,counts[k]))
        #raise Exception()
    ax.legend(plots, labels)
    ax.set_title('factor-target correlations for {0}'.\
                     format(expr_type))

    figtitle = 'ccofs_per_net_{0}_expr'.format(expr_type)
    fpath = figtemplate.format(figtitle)
    f.savefig(fpath)



def make_edge_comparisons(cgraphs, bgraphs):
    cgsets = dict([(k, set(v.edges())) for k, v in cgraphs.iteritems()])

    for bname, bg in bgraphs.iteritems():
        #if bname != 'kn': continue
        f = myplots.fignum(3,(8,8))
        f.clear()
        axes = [f.add_subplot(311),
                f.add_subplot(312),
                f.add_subplot(313)]
        ccolors = dict(zip(cgraphs.keys(), mycolors.getct(len(cgraphs))))

        bgset = set(bg.edges())

        yvals = {'jaccard':[], 'spec':[], 'sens':[]}
        xofs = 0
        heights, xvals ,colors ,names= [], [], [], []
        for cname, cg in sorted(cgraphs.iteritems(),
                                key = lambda x: x[0]):
            cgset = set(cg)
            #SIMILARITIES MATCHING THE ORDER OF SIMNAMES
            yvals['jaccard'].append(float(len(bgset.intersection(cgsets[cname])))/\
                        len(bgset.union(cgsets[cname])))
            yvals['spec'].append(
                    float(len(bgset.intersection(cgsets[cname])))/\
                        len(cgsets[cname]))
            yvals['sens'].append(
                    float(len(bgset.intersection(cgsets[cname])))/\
                        len(bgset))

            #colors.extend([ccolors[cname]] * len(sims))
            #heights.extend(sims)
            names.append(cname )
            #xvals.extend(xofs +arange(len(sims)))
            #xofs = max(xvals) + 2
            #if cname == 'unsup': raise Exception()

        for j, elt in enumerate(yvals.iteritems()):
            metric_name = elt[0]
            heights = elt[1]
            print heights
            ax = axes[j]
            xvals = argsort(argsort(heights))
            ax.bar(xvals, heights, color = [ccolors[n] for n in names])
            ax.set_title('edge similarity vs {0}, metric: {1}'.\
                             format(bname, metric_name))

            myplots.color_legend(f, ccolors.values(), ccolors.keys())
            #for i , n in enumerate(names):
            #    ax.annotate(n, [xvals[i], .001],
            #            xycoords = 'data', 
            #            xytext = [2,0],
            #            textcoords = 'offset points',
            #            rotation = 90, va = 'bottom', ha = 'left')

 
        f.savefig(figtemplate.format('edges_vs_{0}'.format(bname)))
                    

        

def make_transitivity(graphs):
    f = myplots.fignum(3,(8,8))
    udgraphs = dict([(k, nx.Graph(v))
                     for k, v in graphs.iteritems()])
    clusters = dict([(k, nx.algorithms.transitivity(v))
                     for k, v in udgraphs.iteritems()])

    ax = f.add_subplot(111)
    xax = range(len(clusters.keys()))
    ax.plot(xax, clusters.values())
    ax.set_title('transitivity for network graphs')
    ax.set_xticks(xax)
    ax.set_xticklabels(clusters.keys())


    figtitle = 'transitivity'
    fpath = figtemplate.format(figtitle)
    f.savefig(fpath)



def make_clustering(graphs):
    f = myplots.fignum(3,(8,8))

    udgraphs = dict([(k, nx.Graph(v))
                     for k, v in graphs.iteritems()])
    clusters = dict([(k, nx.algorithms.average_clustering(v))
                     for k, v in udgraphs.iteritems()])

    ax = f.add_subplot(111)
    xax = range(len(clusters.keys()))
    ax.plot(xax, clusters.values())
    ax.set_title('clustering for undirected network graphs')
    ax.set_xticks(xax)
    ax.set_xticklabels(clusters.keys())


    figtitle = 'clustering'
    fpath = figtemplate.format(figtitle)
    f.savefig(fpath)



def graph_matrices(seed_graphs, comp_graphs):
    

    for k, v in graphs.iteritems():
        figtitle = '{0}_matrixview'.format()
        fpath = figtemplate.format(figtitle)
        f.savefig(fpath)
