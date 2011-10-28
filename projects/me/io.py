import inference
import networkx as nx 
import cb.utils.memo as mem
def getNet(name = 'easy0', **kwargs):
    if name == 'easy0':
        return inference.get_easy0(**kwargs)
    else: raise Exception()
    return 



def getGraph(name = 'easy0', control = 'in_degree',atype='wormtile',**kwargs):
    edges = getNet(**mem.rc(kwargs, atype = atype))
    dg = nx.DiGraph()
    dg.add_weighted_edges_from(edges)

    if control == 'in_degree':
        in_degrees = dict([(n,0) for n in dg.nodes()])
        for e in dg.edges(): in_degrees[e[1]] += 1
        for n in dg.nodes():
            if len(dg[n]) == 0 and in_degrees[n] < 4:
                dg.remove_node(n)
    
    return dg
