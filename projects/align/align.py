import cb.utils.memo as mem
import cb.p.network.io as nio
import cb.p.nfilt.comparator as comp
import cb.p.brainmap.worm as worm
import cb.config as cfg

import cb.utils.graphs.utils as gu
import cb.utils.graphs.draw as gd
import cb.utils.plots as myplots

import matplotlib.pyplot as plt
from numpy import *
import  numpy as np
import os, re, itertools as it
import subprocess as spc

import networkx as nx

def fly_graphs(params):
    '''
params:

sort_criterion: 'total_degree', 'out_degree'
graphs:         'b_v_m', 'b_v_b'
edge_frac:      ...
node_frac:      ...
max_nodes:      200
'''
    all_graphs = comp.get_graphs(restriction = 'none')   
    gkeys = params.get('nets',['bn','mn'])
    bnet = all_graphs[gkeys[0]]
    mnet = all_graphs[gkeys[1]]
    g0s = [bnet,mnet]

    max_nodes = params.get('max_nodes', 50)
    sort_criterion = params.get('sort', 'out_degree')
    
    nodes = set(mnet.nodes()).intersection(bnet.nodes())
    g0s = [gu.restricted_graph(g,list(nodes)) for g in g0s]

    outs = [ dict([(n,len(g[n])) for n in nodes])
             for g in g0s]
    ins = [dict([(n,0) for n in nodes]) for g in g0s]
    for i,g in enumerate(g0s):
        for e in g.edges(): 
            ins[i][e[1]] += 1

    #if sort_criterion == 'total_degree':
    #     srtnodes = sorted(nodes,key =  lambda x: min(outs[ + ins[x])[::-1]
    if sort_criterion == 'out_degree':
        srtnodes = sorted(nodes, 
                          key = lambda x: min(outs[0][x],outs[1][x]) )[::-1]

    nodelist = srtnodes[:max_nodes]
    graphs = []
    for g0 in [all_graphs[gkeys[0]],all_graphs[gkeys[1]]]:
        g0 = gu.restricted_graph(g0, nodelist)
        #compute a noisy verson of the binding network
        ge = g0.edges()
        rnd_edges = set([ge[i] for i in random.permutation(len(ge))\
                             [:len(ge) * params.get('edge_frac',1)]])    
        n0 = g0.nodes()
        #choose a random sampling of node according to node frac
        rnd_nodes = [n0[i] for i in random.permutation(len(g0))\
                         [:len(g0) * params.get('node_frac',1)]]
        #dummy run with no parameters

        g1 = gu.restricted_graph(g0, rnd_nodes)
        for e in g1.edges(): 
            if not e in rnd_edges: g1.remove_edge(*e)
        graphs.append(g1)
        


    return graphs
    
def brain_graphs(params):
    g0 = gu.restricted_graph(worm.get_graph(), 
                             worm.get_graph().nodes()[:])
    n0 = g0.nodes()
    #choose a random sampling of node according to node frac
    rnd_nodes = [n0[i] for i in random.permutation(len(g0))\
                     [:len(g0) * params['node_frac']]]
    ge = g0.edges()
    rnd_edges = set([ge[i] for i in random.permutation(len(ge))\
                     [:len(ge) * params['edge_frac']]])
    
    #dummy run with no parameters
    g1 = gu.restricted_graph(g0, rnd_nodes)
    for e in g1.edges(): 
        if not e in rnd_edges: g1.remove_edge(*e)

    graphs = [g0,g1]
    return graphs

def fly_params_list():
    nets = ['kn','mn','bn']
    params = []
    for net0 in nets:
        for net1 in nets:
            for edge_frac in [.9]:
                for max_nodes in [50]:
                    params.append(
                        dict(nets = [net0, net1],
                             edge_frac = edge_frac,
                             graphs = 'b_v_m',
                             max_nodes = max_nodes,
                             name = '{0}_v_{1};_efrac_{2};{3}_nodes'\
                                 .format(net0,net1,edge_frac,max_nodes)))
    return params
                    

def get_params(name, p):
    if name == 'fly':
        return fly_params_list()[p]
    elif name == 'brain':
        return [{'node_frac':1,
                 'edge_frac':1,
                 'name':'identity'},
                {'node_frac':.9,
                 'edge_frac':1,
                 'name':'nodes_missing'},
                {'node_frac':1,
                 'edge_frac':.9,
                 'name':'edges_missing'},
                {'node_frac':.75,
                 'edge_frac':1,
                 'name':'more_nodes_missing'},
                {'node_frac':1,
                 'edge_frac':.75,
                 'name':'more_edges_missing'}][p]
        

def align_graphm(name = 'brain',
                 do_plot_initial = True,
                 p = 0):
    params = get_params(name, p)

    if name == 'fly':
        graphs = fly_graphs(params)
    elif name == 'brain':
        graphs = brain_graphs(params)

    graphs = [nx.convert.convert_to_undirected(i) for i in graphs]
    
    
    g0_tmp = nx.Graph(graphs[0])
    for n in graphs[1].nodes(): 
        if not g0_tmp.has_node(n): g0_tmp.add_node(n)

    embeddings =[gd.getpos(g) for g in [g0_tmp, graphs[1]]]
    if do_plot_initial:
                     plot_initial(graphs,name,embeddings,
                                  pname = params.get('name', 'generic'))

    root = cfg.dataPath('align/run_{0}_pset{1}'.format(name,p))
    if not os.path.isdir(root): os.makedirs(root)
    gfiles = [os.path.join(root, 'graph_0'),
              os.path.join(root, 'graph_1')]
    
    array_keys = [sorted(g.nodes()) for g in graphs]
    gmats = [array(nx.convert.to_numpy_matrix(g, nodelist = array_keys[i]))
             for i, g in enumerate(graphs)]
    #labels = [nx.convert.convert_node_labels_to_integers(g)
    #          for g in graphs]

    gascii= ['\n'.join([' '.join(['{0}'.format(1 if e > 0 else 0) 
                                   for e in row]) 
                        for row in graph]) + '\n'
             for graph in gmats]

    for i, ga in enumerate(gascii):
        f = open(gfiles[i], 'w')
        f.write(ga)
        f.close()
              
    sim_file = os.path.join(root, 'sims_01')
    simtype = 'null'
    if simtype != 'null':
        sims = compute_sims(graphs, 'null')
        sf = open(sim_file, 'w')
        sf.write( '\n'.join([' '.join(['{0}'.format(e) 
                                       for e in row]) 
                             for row in sims]))
    else:
        if os.path.isfile(sim_file):
               os.remove(sim_file)
    outfile = os.path.join(root, 'out')
    if os.path.isfile(outfile): os.remove(outfile)
    ca = config_ascii( gfiles, sim_file, outfile)
    configfile = os.path.join(root, 'config.txt')
    cf = open(configfile, 'w')
    cf.write(ca)
    cf.close()

    cmd = 'graphm {0}'.format(configfile)
    print cmd
    prc = spc.Popen(cmd,
                    shell = True,
                    stdout = spc.PIPE)
    out = prc.stdout.read()
    
    output = {'graphs': graphs,
              'outfile': outfile,
              'run_output':out,
              'cmd':cmd,
              'name':name,
              'params':params,
              'keys': array_keys,
              'embeddings':embeddings}
    import pickle
    pickle.dump(output, open(os.path.join(root, 'data.pickle'), 'w'))
    

    return output

def view(output):
    outfile = output['outfile']
    graphs = output['graphs']

    results = open(outfile).readlines()
    
    read_cols = 0
    read_perms = 0
    perms = []
    for r in results:
        if r[0:12] == 'Permutations':
            read_cols = 1
            continue
        elif  read_cols:
            cols = re.compile('\s+').split(r.strip())
            read_cols = 0; read_perms =1;
            continue
        elif read_perms:
            perms.append(re.compile('\s+').split(r.strip()))
    perms = array(perms, int) - 1
    
    matched = plot_perms(perms,output)
    return perms  , matched

def config_ascii(gfiles, simfile, outfile):
    return '''//*********************GRAPHS**********************************
//graph_1,graph_2 are graph adjacency matrices, 
//C_matrix is the matrix of local similarities  between vertices of graph_1 and graph_2. 
//If graph_1 is NxN and graph_2 is MxM then C_matrix should be NxM
graph_1={0} s
graph_2={1} s
C_matrix={2} s
//*******************ALGORITHMS********************************
//used algorithms and what should be used as initial solution in corresponding algorithms
algo=I U RANK QCV rand PATH s
algo_init_sol=unif unif unif unif unif unif s
solution_file=solution_im.txt s
//coeficient of linear combination between (1-alpha_ldh)*||graph_1-P*graph_2*P^T||^2_F +alpha_ldh*C_matrix 
alpha_ldh=0 d
cdesc_matrix=A c
cscore_matrix=A c
//**************PARAMETERS SECTION*****************************
hungarian_max=10000 d
algo_fw_xeps=0.01 d
algo_fw_feps=0.01 d
//0 - just add a set of isolated nodes to the smallest graph, 1 - double size 
dummy_nodes=0 i
// fill for dummy nodes (0.5 - these nodes will be connected with all other by edges of weight 0.5(min_weight+max_weight))
dummy_nodes_fill=0 d
// fill for linear matrix C, usually that's the minimum (dummy_nodes_c_coef=0),
// but may be the maximum (dummy_nodes_c_coef=1)
dummy_nodes_c_coef=0.01 d

qcvqcc_lambda_M=10 d
qcvqcc_lambda_min=1e-5 d


//0 - all matching are possible, 1-only matching with positive local similarity are possible
blast_match=1 i
blast_match_proj=0 i


//****************OUTPUT***************************************
//output file and its format 
exp_out_file={3} s
exp_out_format=Parameters Compact Permutation s
//other
debugprint=0 i
debugprint_file=debug.txt s
verbose_mode=1 i
//verbose file may be a file or just a screen:cout
verbose_file=cout s
'''.format(gfiles[0], gfiles[1], simfile, outfile)

def compute_sims(graphs, sim_type):
    if sim_type == 'null':
        return ones((len(graphs[0]), len(graphs[1])))
    else:
        raise Exception('sim_type {0} not yet implemented'.\
                            format(sim_type))

def plot_perms(perms, out):
    embeddings = out['embeddings']
    keys = out['keys']
    name = out['name']
    pname =out['params']['name'] 

    graphs = out['graphs']
    all_nodes =set( list(it.chain(*out['keys'])))
    p0 = [p[0] for p in perms]
    p5 = [p[3] for p in perms]

    nodes = out['keys']
    matching = [(nodes[0][p0[i]] if p0[i] < len(nodes[0]) else -1, 
                 [
                nodes[1][perms[i][j]]
                if perms[i][j] < len(nodes[1]) else -1
                for j in range(1,6)
                ])
                for i in range(len(p0))]
                

    f = myplots.fignum(3, (8,5))
    
    bpos, mpos = embeddings
    
    all_ax = [121, 222, 224]
    all_graphs = [graphs[0],graphs[1],graphs[1]]
    all_pos = [bpos,mpos,bpos]
    all_names = ['{0} net 0; red=qp, orange=other'.\
                      format(name),
                 '{0} net 1, type is {1}'.\
                      format(name, pname),
                 '{0} net 1, net0 projection'.\
                      format(name, pname)]
                 
    alpha = .05
    n_any_matches = dict([(e[0], True if e[0] in e[1] else False)
                          for e in matching])

    n_4_matches = dict([(e[0], True if e[0] == e[1][2] else False)
                          for e in matching])


    for a,g,p,n in zip(all_ax,all_graphs,all_pos,all_names):
        ax = f.add_subplot(a)
        ax.set_title(n)
        plt.gca().set_xticks([]); plt.gca().set_yticks([])
        gd.easy_draw(g,p, alpha = alpha)
    

        nodes_4_matched =  [n for n in g.nodes() if n_4_matches.get(n,False)]
        nodes_matched =  [n for n in g.nodes() 
                          if n_any_matches.get(n,False) 
                          and not n in nodes_4_matched]

        gd.draw(g, p,g.edges(),
                skw = dict(color = 'orange',
                           s = 20,
                           alpha = 1),
                scatter_nodes = nodes_matched)

        gd.draw(g, p,g.edges(),
                skw = dict(color = 'red',
                           s = 20),
                scatter_nodes = nodes_4_matched)


    f.savefig(myplots.figpath('matched_graphs_{0}_p{1}.pdf'.format(name,pname)))

    return matching
    
def plot_initial(graphs, name, embeddings, pname = 'generic'):
    f = myplots.fignum(3, (8,4))
    
    bpos, mpos = embeddings

    ax0 = f.add_subplot(121)
    plt.gca().set_xticks([]); plt.gca().set_yticks([])
    ax0.set_title('{0} network 0'.\
                      format(name))
    gd.easy_draw(graphs[0], bpos)

    ax1 = f.add_subplot(222)
    ax1.set_title('{0} network 1, type is {1}'.\
                      format(name, pname))
    plt.gca().set_xticks([]); plt.gca().set_yticks([])
    gd.easy_draw(graphs[1], mpos)


    ax2 = f.add_subplot(224)
    ax2.set_title('{0} net 1, net0 projection'.\
                      format(name, pname))
    gd.easy_draw(graphs[1], bpos)
    plt.gca().set_xticks([]); plt.gca().set_yticks([])

    f.savefig(myplots.figpath('basic_graphs_{0}_p{1}.pdf'.format(name,pname)))

    return [bpos, mpos]
