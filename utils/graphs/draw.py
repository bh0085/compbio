import matplotlib.pyplot as plt
from numpy import *
import networkx as nx
from matplotlib.pyplot import gca

def getpos(g):
    return nx.graphviz_layout(g)

def easy_draw(graph,
              pos,
              alpha = .25):
    draw(graph, pos, graph.edges(),
         skw = {'s':10,
                'edgecolor':'black',
                'facecolor':'white'},
         
         ckw = dict([(k,dict(color = 'black',
                             alpha = alpha,
                             linewidth = 1,
                             arrowstyle = '->'))
                     for k in graph.edges()]))
    gca().set_xticks([])
    gca().set_yticks([])

def draw(graph, pos,edges,
         labels = None,
         scatter_nodes = None,
         ckw = {},     #CONNECTION STYLE KWS
         skw = {}):    #SKATTERPLOT KWS):
    '''
Draw a graph with scatterpoint keyword given in skw
and connectionpoint keyworks given by ckw.

ckw:
  ckw takes the form of a dictionary; edges lacking a key in ckw
  will not be plotted. For every edges that is plotted, the following
  keywords are allowed:

  'alpha'
  'fc'
  'ec'

  skw:
  see pyplot.scatter?

'''
    #MATPLOTLIB ARROW TYPES
    '''
-	None
->	head_length=0.4,head_width=0.2
-[	widthB=1.0,lengthB=0.2,angleB=None
-|>	head_length=0.4,head_width=0.2
<-	head_length=0.4,head_width=0.2
<->	head_length=0.4,head_width=0.2
<|-	head_length=0.4,head_width=0.2
<|-|>	head_length=0.4,head_width=0.2
fancy	head_length=0.4,head_width=0.4,tail_width=0.4
simple	head_length=0.5,head_width=0.5,tail_width=0.2
wedge	tail_width=0.3,shrink_factor=0.5
'''
    #cstr = "arc3,rad="+'.2'#str(-(i - nt/2) *.2)

    if scatter_nodes == None:
        scatter_nodes = graph.nodes()

    ax = plt.gca()

    for i, e in enumerate(edges):
        if e in ckw:

            
            ax.annotate('', xy=pos[e[1]],  xycoords='data',
                        xytext=pos[e[0]], textcoords='data',
                        arrowprops=dict( 
                                         #connectionstyle=cstr,
                                         #shrinkA=10,shrinkB = 10,
                                         **ckw.get(e,{})  ),
                        )

    
    #return
    xys  = [pos[n] for n in scatter_nodes]
    if xys:
        ax.scatter(*zip(*xys), **skw)
    
def overlay(graph, pos, edges, 
            colors = {},alphas = {}, tripoints = {}, 
            tristretch = .5, #TRIPOINTS KW
            circlepoints = {}, circleradii = {}    #CIRCPOINTS KW
            ):
    '''
Overlay data on a subset of graph edges.

A few modes are so far supported:
1) Tripoints:
   Stretches a triangle off of each edge towards the 
   directing point given as a 2vector in tripoints.

2) Circlepoints:
   A list of dicts keyed (optionally) with
     lw: linewidth
     ec: edgecolor
     s:  size
     
   And colored by the (shared) colors key
'''
    ax = plt.gca()
    if colors == None:
        colors = zeros((len(edges),3))

    for i, e in enumerate(edges):
        xy0 = array(pos[e[0]])
        xy1 = array(pos[e[1]])
        if tripoints.has_key(e):
            xy2 = (xy0 + xy1)/2 +\
                (tripoints.get(e) - (xy0 + xy1)/2) * tristretch
            from matplotlib.patches import Polygon
            poly = Polygon([xy0,xy1,xy2], 
                           color = colors.get(e,'black'),
                           alpha = alphas.get(e,1))
            ax.add_artist(poly)
            
        
        if circlepoints.has_key(e):
            xyc =(xy0+xy1)/2
            cp = circlepoints[e]

            import matplotlib.transforms as transforms
            transform = transforms.ScaledTranslation(\
                xyc[0],xyc[1],ax.transData )

            from matplotlib.patches import Circle
            poly = Circle([0,0], circleradii.get(e,10),
                          transform = transform,
                          **cp)
            ax.add_artist(poly)

            #poly = Circle([5,5], circleradii.get(e,10),
            #              color = 'black',
            #              transform = transform,
            #              alpha = .1)#**cp)
            #ax.add_artist(poly)
            #ckw.append(cp)
            #colc.append( colors[i])
            #lwc.append(cp.get('lw',1))
            #ecc.append(cp.get('ec',[0,0,0]))
            #sc.append(cp.get('s', 100))
            #linestyles
             
