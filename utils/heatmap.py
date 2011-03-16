'''heatmap v. .1. 
 Call heatmap.heatmap() to get a heatmap.
'''

import compbio.utils.colors as mycolors
import matplotlib.pyplot as plt
from numpy import *
import numpy as np
def heatMap(grid, 
            ann = [],
            xlabel = 'none', 
            ylabel = 'none',
            **kwargs):

  '''
  for now, this is a primordial heatmap script.
  The optimal use case is as in the function heatmapGene from 
  compbio.projects.predict.py. In this case, we input a grid and
  a list of annotations that happen to annotate each and every
  grid point. In particular, the dictionary contains an xvalue 
  and a yvalue for each point as well as an entry, 'pkeys' that
  give the names of the elements being plotted on the x and the y.

  e.g: ann[0] = {'pvalue':0.,
                 'cvalue':0.,
                 'pkeys':['pvalue','cvalue']
                 }
                 
  Of course, this assumes a standard form for annotations.
  I guess its not so bad for now.

'''

  cmap = mycolors.blackbody()

  ax = kwargs.get('axes', None)
  if not ax:
    f = kwargs.get('fig', 0)
    figure = plt.figure(f)
    ax = figure.add_subplot(111)
    
  ax.imshow(grid.T, cmap = cmap, interpolation = 'nearest',
            origin = 'lower', aspect = 'auto')

  xticks, xticklabels, yticklabels = [], [], []
  yticks = []

  edgept = lambda i,j:  ( i == 0 and j == 0) \
      or ( i == len(ann) -1 and j == len(ann[0])-1)
  lowpt = argmin(grid)
  highpt = argmax(grid)

  ishigh = lambda i,j : (i,j) == unravel_index(highpt, shape(grid))
  islow = lambda i,j: (i,j) == unravel_index(lowpt, shape(grid))


  
  for i in range(shape(ann)[0]):
    for j in range(shape(ann)[1]): 
      
      if random.random() < ( 0./product(shape(ann))) \
            or edgept(i,j) or islow(i,j) or ishigh(i,j): 

        d = ann[i][j]
        dkeys = d['pkeys']
        s =  dictString(ann[i][j])
        xy = [i,j]
        color = 'black'

        xytext=(30,10)
        textcoords='offset pixels'
        bbox = None #bbox=dict(boxstyle="round4", fc="none")
        arrowprops=dict(arrowstyle="-|>",
                        connectionstyle="arc3,rad=-0.2",
                        relpos=(0., 0.),
                        shrinkA = 0,
                        shrinkB = 10,
                        fc="none") 
        ax.scatter(xy[0],xy[1], 100, 
                    color = 'none',
                    edgecolor = 'black')
        if islow(i,j): 
          color ='red'
          s+= '\nLow Point: {0}'.format(strFormat(grid[i,j]))
        elif ishigh(i,j):
          color = 'blue'
          s+= '\nHigh Point: {0}'.format(strFormat(grid[i,j]))


        if xlabel == 'none': xlabel = dkeys[0]
        if ylabel == 'none': ylabel = dkeys[1]

        if not edgept(i,j): ax.annotate(s,xy,xytext = xytext,
                                        textcoords = textcoords,
                                        bbox = bbox,
                                        arrowprops =arrowprops)
        
        xticks.append(i)
        xticklabels.append(strFormat(d[dkeys[0]]))
        yticks.append(j)
        yticklabels.append(strFormat(d[dkeys[1]]))
        
  ax.set_xticks(xticks)
  ax.set_xticklabels(xticklabels)
  ax.set_yticks(yticks)
  ax.set_yticklabels(yticklabels)
  ax.set_xlabel(xlabel)
  ax.set_ylabel(ylabel)

  return ax

def dictString(d):
  #s0 = ''
  #for k in d['pkeys']:
  #  v = d[k]
  #  s0 += '{0}: {1}\n'.format(k, strFormat(v)) % v
  s0 = '('+', '.join([strFormat(d[k]) for k in d['pkeys']]) + ')'
  return s0

def strFormat(f, e_range = [-3,3]):
  
  
  if log10(f) > e_range[0] and log10(f) < e_range[1]:
    return '%3.2f' % f
  else:
    return '%1.2e' % f
  
