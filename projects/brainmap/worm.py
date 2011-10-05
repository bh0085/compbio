import os, itertools as it
import cb.config as cfg
import cb.utils.plots as myplots
import cb.utils.colors as mycolors
import xlrd 
import networkx as nx
from numpy import *
import numpy as np
import cb.utils.graphs.draw as gd
import cb.utils.memo as mem

'''NeuronConnect:

Type: Type of synapse: S: Send or output (Neuron 1 pre-synaptic to Neuron 2); 

Sp: Send-poly (Neuron 1 is pre-synaptic to more than one postsynaptic partner.  
Neuron 2 is just one of these post-synaptic neurons, see Figure 1 below.  
In White et al, 1986, these polyadic synaptic connections were denoted by 
'm' in the tables of Appendix 1); 

R: Receive or input (Neuron 1 is post-synaptic to Neuron 2); 

Rp: Receive-poly (Neuron 1 is one of several post-synaptic partners of Neuron 2.  
See Figure 1 and above); 

EJ: Electric junction; 
NMJ: Neuromuscular junction (only reconstructed NMJ's are represented). 
'''
defplots = {}
def load(plots = defplots,
         reset = False):
    kwargs = dict(reset = reset)
    edge_set = get_edge_set()
    g = get_graph(**mem.sr(kwargs))
    pos = get_pos(**mem.sr(kwargs))
    trips = set([])
    for k1 in g:
        for k2 in g[k1].keys():
            for k3 in g[k1].keys():
                if g[k2].has_key(k3):
                    trips.add((k2,k3,k1))

    tripoints = dict([((e[0],e[1]),pos[e[2]]) for e in trips])
    if plots.get('basic_structure', False):
        f = myplots.fignum(1)
        ax = f.add_subplot(111)
        gd.draw(g, pos, g.edges(),
                skw = {'s':10,
                       'edgecolor':'black',
                       'facecolor':'white'},

                ckw = dict([(k,dict(color = 'black',
                                    alpha = .25,
                                    linewidth = 1,
                                    arrowstyle = '->'))
                            for k in g.edges()]))
        
        f.savefig(myplots.figpath('basic_structure_edges={0}'.format(edge_set))) 

    if plots.get('feed_forward', True):
        gd.overlay(g,pos,g.edges(),
                   tripoints = tripoints, 
                   alphas = dict([(e,.1) for e in g.edges()]))
                   
        f.savefig(myplots.figpath('feed_forward_edges={0}'.format(edge_set)))

    if plots.get('degrees' , False):
        make_degree_plots_0();

    maxflow = nx.algorithms.ford_fulkerson(g, 'AVAL','PVPL','weight')

    imaps = get_array_imaps()
    nnames = imaps['nnames']
    node_stats = dict([(k,{}) for k in nnames])
    for k,v in node_stats.iteritems():
        v['out_degree'] = len([e for e in g.edges() if e[0] == k])
        v['in_degree'] = len([e for e in g.edges() if e[1] == k])
        
    f = myplots.fignum(3, (12,6))
    outs = [v['out_degree'] for k, v in node_stats.iteritems()]
    ins =[v['in_degree'] for k , v in node_stats.iteritems()]
    raw_data= array([outs,ins]).T





    make_data_transform(raw_data)
    data = transform_data(raw_data)
    kd = make_kdtree(data)
    k = 5 
    nn = compute_nns(kd, k)
    knn= nn['nn']
    knn_dists = nn['dists']
    dists = compute_dists(data)


    mean_dists = np.mean(knn_dists[:,1:],1)
    mean_colors =sqrt(mean_dists[:,newaxis] * [1/np.max(mean_dists), 0,0])

    ax = f.add_subplot(121)

    ax.scatter(data[:,0],data[:,1],s = 15,
              facecolor = mean_colors,
              edgecolor = 'none'
               )
    
    ax.set_xlabel('scaled out degree')
    ax.set_ylabel('scaled in degree')

    ax2 = f.add_subplot(122)

    ax2.imshow(dists,
               interpolation = 'nearest',
               aspect = 'auto')
    ax2.set_title('distance matrix for scaled degrees')
    
    f.savefig(myplots.figpath('distances_{0}'.format(edge_set)))
        

    return g

data_means = None
data_scales= None
def make_kdtree(data):
    from scipy.spatial import KDTree
    kd = KDTree(data)
    return kd 

def compute_nns(kd,k):
    query = kd.query(kd.data,k+1)
    return {'nn':query[1],
            'dists':query[0]}

def make_data_transform(data):
    globals()['data_means'] = mean(data,0)
    globals()['data_scales'] = var(data - data_means[newaxis,:], 0)
    
def transform_data(data):
    return ( data - data_means[newaxis,:] ) / data_scales[newaxis,:]

def compute_dists(data):
    dotprod = np.sum(data[:,newaxis,:]*data[newaxis,:,:],2)
    n = len(dotprod)
    sqrs = resize(dotprod.diagonal(), (n,n))
    dists = sqrs + sqrs.T - 2 * dotprod
    return dists

def get_edge_set():
    return ('S',)
def get_rows(**kwargs):
    def set_rows(**kwargs):
        root = cfg.dataPath('wormbrain/2006')
        connect_file = os.path.join(root, 'NeuronConnect.xls')
        fp_file = os.path.join(root,'NeuronFixedPoints.xls')
       
        cwb = xlrd.open_workbook(connect_file)
        sh = cwb.sheets()[0]
       
        rows = [[e.value for e in sh.row(i)] for i in range(1,sh.nrows) ]
        return rows
    return mem.getOrSet(set_rows, **mem.rc(kwargs))
        
def get_synapse_dict(**kwargs):
    def set_synapse_dict(**kwargs):
        rows = get_rows()
        all_out_cxns =  dict([(k, [ e for e in list(val) ])
                           for k,val in it.groupby(\
                   sorted(rows, key = lambda x: x[0]),
                   key = lambda x: x[0])])
        return all_out_cxns
    return mem.getOrSet(set_synapse_dict,
                        **mem.rc(kwargs))
       

def get_synapse_array(**kwargs):
    def set_synapse_array(**kwargs):

        imaps = get_array_imaps(**mem.sr(kwargs))
        ctypes =imaps['ctypes']
        ctypes_imap = imaps['ctypes_imap']
        nnames = imaps['nnames']
        nnames_imap = imaps['nnames_imap']
        
        all_out_cxns = get_synapse_dict(**mem.sr(kwargs))

        cxns = zeros((len(nnames),len(nnames), len(ctypes)))
        for k1,rows in all_out_cxns.iteritems():
           for row in rows:
               cxns[nnames_imap[k1],nnames_imap[row[1]],
                    ctypes_imap[row[2]]] +=row[3]  

        return cxns
    return mem.getOrSet(set_synapse_array,
                        **mem.rc(kwargs))
def get_array_imaps(**kwargs):
    def set_array_imaps(**kwargs):
        sdict =get_synapse_dict(**mem.sr(kwargs))
        nameset = set([])
        for k,v in sdict.iteritems():
            nameset.add(k)
            nameset.update([r[1] for r in v])
        nnames = list(nameset)

        ctypes = [u'Rp', u'EJ', u'Sp', u'S', u'R', u'NMJ']
        ctypes_imap = dict([(k,i) for i, k in enumerate(ctypes)])
        nnames_imap = dict([(k,i) for i, k in enumerate(nnames)])

        return {'ctypes':ctypes,
                'ctypes_imap':ctypes_imap,
                'nnames':nnames,
                'nnames_imap':nnames_imap}
    return mem.getOrSet(set_array_imaps,
                        **mem.rc(kwargs))


def get_pos(**kwargs):
    def set_pos(**kwargs):
        g = get_graph(**mem.sr(kwargs))
        pos = gd.getpos(g)
        return pos
    return mem.getOrSet(set_pos, **mem.rc(kwargs))

def get_graph(**kwargs):
    def set_graph(**kwargs):
        edge_set = get_edge_set()
        rows = get_rows(**mem.sr(kwargs))
        sub_cxns = dict([(k, [ e for e in list(val) if e[2] in edge_set])
                         for k,val in it.groupby(\
                    sorted(rows, key = lambda x: x[0]),
                    key = lambda x: x[0])])
        g = nx.DiGraph();
        for k, v in  sub_cxns.iteritems():
            g.add_weighted_edges_from([(e[0], e[1], e[3]) for e in v])
        return g
    return mem.getOrSet(set_graph,**mem.rc(kwargs))
        


def make_degree_plots_0():
        cxns = get_synapse_array()
        rows = get_rows()

        imaps = get_array_imaps()
        ctypes =imaps['ctypes']
        ctypes_imap = imaps['ctypes_imap']
        nnames = imaps['nnames']
        nnames_imap = imaps['nnames_imap']
            
        f2 = myplots.fignum(2, (12,6))
        ax1 = f2.add_subplot(121)
        ax2 = f2.add_subplot(122)
        
        var_degs = np.sum(cxns,1)
        maxval = log10(np.max(var_degs) + 1)
        ct = mycolors.getct(len(ctypes))
        for z in range(len(ctypes)):
            vals = var_degs[:,z]
            vals = log10(1 + vals)
            count,kde = make_kde(vals)
            xax = linspace(0,maxval,10)
            h = histogram(vals, xax)
            ax1.hist(vals,xax, 
                     color = ct[z],
                     zorder = 10,
                     alpha = .25)
            ax1.plot(xax,kde(xax)*sum(h[0]),
                             label = ctypes[z],
                             color = ct[z],
                             zorder = 5)
            ax1.set_xlabel('$log_10$ of edge degrees of various types')
        ax1.legend()
        
        logxy = [ log10(1 +var_degs[:,ctypes_imap['S']]),
                  log10(1 +var_degs[:,ctypes_imap['R']])]
        max_inode =np.argmax(logxy[0] + logxy[1])
        max_nodename = [k 
                        for k,v in nnames_imap.iteritems() 
                        if v == max_inode][0]
        

        ax2.scatter(logxy[0]+.15*random.rand(len(nnames))
                    ,logxy[1] + .15*random.rand(len(nnames)),
                    color = 'red',
                    alpha = .3)
        ax2.set_xlabel('Sending Degree')
        ax2.set_ylabel('Receiving Degree')
        r2 = corrcoef(logxy[0],logxy[1])[1,0]

        myplots.maketitle(ax2, ('correlation coeff: {0:2.2},\n'+\
                              'max {1} has {2} $e_{{out}}$, {3} $e_{{in}}$')\
                              .format(r2, max_nodename, 
                                      var_degs[max_inode, ctypes_imap['S']],
                                      var_degs[max_inode, ctypes_imap['R']]))
        myplots.maketitle(ax1, 'histogram and KDE of\nvarious edge degrees')        
        f2.savefig(myplots.figpath('degree_histograms_{0}'.format(edge_set)))


def make_kde(vals):
    from scipy.stats import gaussian_kde
    count = len(vals)
    dist = array(vals)
    density = gaussian_kde(dist)
    density.covariance_factor = lambda : .5
    density._compute_covariance()
    kde = (density)
    
    return count, kde
