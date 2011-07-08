import cb.utils.plots as myplots
import cb.utils.colors as mycolors
import networkx as nx

import matplotlib.pyplot as plt

def show(g, pos = None, **kwargs):
    if pos == None:
        pos = nx.graphviz_layout(g)
    nx.draw(g,pos, 
            arrowstyle="->",
            connectionstyle="angle3,angleA=0,angleB=-90",
            **kwargs)
    return

def show_sims(sims):
    ax = plt.gca()
    ax.imshow(sims, cmap = mycolors.blackbody(),interpolation = 'nearest')
