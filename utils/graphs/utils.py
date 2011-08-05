import networkx as nx
import cb.utils.memo as mem

def consensus_graph(name = 'none',
                    graphs =(),
                    **kwargs):
  '''
Get a graph having nodes consisting of the union o
the nodes in all graphs and having edges consisting
of the intersection of edges in all graphs.

'''
  def get_cons_graph(**kwargs):
    graphs = kwargs.get('graphs')
    if type(graphs[0]) != nx.DiGraph:
      raise Exception('For now, this method is only compatible with digraph')
    all_nodes = set.union(*[set(g.nodes()) for g in graphs])
    ##NOTE, THIS SYNTAX IS DESIGNED FOR DIRECTED GRAPHS
    ##FOR UNDIRECTED, IT WILL FAIL TO COUNT BIDIRECTIONAL EDGES
    all_edges = set.intersection(*[set(g.edges()) for g in graphs])
    cons = nx.DiGraph()
    cons.add_nodes_from(all_nodes)
    cons.add_edges_from(all_edges)

    return cons
  return mem.getOrSet( get_cons_graph, 
                      **mem.rc(kwargs,
                               on_fail = 'compute',
                               register =name ,
                               graphs = graphs,
                               name = name))

def restricted_graph(graph, nodelist):
    gnew = type(graph)(graph)
    diff_nodes = set.difference( set(graph.nodes()), set(nodelist))
    gnew.remove_nodes_from(diff_nodes)
    
    return gnew
