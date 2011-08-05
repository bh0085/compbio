import cb.p.network.utils as nu
import cb.p.nfilt.comparator as comp
import cb.p.network.io as nio
import cb.utils.graphs.utils as gu
import cb.utils.graphs.draw as gd
import cb.utils.plots as myplots
import networkx as nx



def run0(nodelist = 'bdtnp',
         term_layout = False,
         demand = 'motif'):
    '''
Starting with a network graph, grab the controlled flybase terms
for each node. If n_terms is specified, use only the top 'N'
recurring terms.

Plot the nodes and their connections on the graph layout from the
controlled vocabulary graph layout.

kwargs:
    node_list:     how to draw the nodes to compare. ('kn', 'bdtnp')

'''
    graphs = comp.get_graphs()

    nterms = 30
    term_groups = nu.term_groups(nodelist, nterms = nterms)
    term_net = nu.term_network(nodelist, nterms = nterms)

    m_b_cons = gu.consensus_graph(name = 'binding_motif',
                                 graphs = (graphs['bn'],graphs['mn'])
                                 )


    m_b_k_cons = gu.consensus_graph(name = 'binding_motif_knowledge',
                                 graphs = (graphs['bn'],graphs['mn'],
                                           graphs['kn'])
                                 )

    
    if demand == 'knowledge':  reg_net =  gu.restricted_graph(m_b_k_cons, nio.getBDTNP().keys())
    elif demand == 'binding': reg_net = gu.restricted_graph(graphs['bn'], nio.getBDTNP().keys())
    elif demand == 'motif': reg_net  = gu.restricted_graph(m_b_cons, nio.getBDTNP().keys())
    else: raise Exception()

    ctd_nodes =list( set.union(set([e[0] for e in term_net.edges()]),
                          set([e[1] for e in term_net.edges()]),
                          set([e[0] for e in reg_net.edges()]),
                          set([e[1] for e in reg_net.edges()])))

    reg_net = gu.restricted_graph(reg_net, ctd_nodes)
    reg_net.add_nodes_from(term_net.nodes())
    nodes = reg_net.nodes()
    edges = reg_net.edges()

    
    

    
    pos = nx.graphviz_layout(reg_net) if not term_layout\
        else nx.graphviz_layout(term_net)
    
    fig = myplots.fignum(3, (8,8))
    ax = fig.add_subplot(111)
    

    semantic_edges = term_net.edges()
    reg_edges = reg_net.edges()
    reg_in = set([e[0] for e in reg_edges])
    reg_out = set([e[1] for e in reg_edges])
    semantic_plotted_edges = [e for e in semantic_edges
                              if e[0] in reg_in or e[1] in reg_in ]
    

    gd.draw(reg_net,pos,reg_edges,
            scatter_nodes = nodes,
            skw = {'color':'none', 'edgecolor':'none'},
            ckw = dict([[e, {'alpha' :.5,
                             'ec' : 'red',
                             'arrowstyle': '->',
                             'zorder': 2}] 
                        for e in reg_edges]))

    gd.draw(reg_net,pos,semantic_plotted_edges, 
            scatter_nodes =[],
            ckw = dict([[e, {'alpha' :.5,
                             'ec' : 'black',
                             'arrowstyle': '|-|',
                             'fc' : 'blue' }] for e in semantic_plotted_edges]))


    
    path = myplots.figpath('{0}.pdf'.\
                               format('{0}_terms={1}{2}'.
                                      format('term_layout' if term_layout
                                             else 'regulation_layout',
                                             nterms,
                                             demand)))    
    fig.savefig(path)

    raise Exception()
