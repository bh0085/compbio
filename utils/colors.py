import random
from numpy import *
import numpy as np
import matplotlib.colors as mplcolors
import matplotlib
import matplotlib.pyplot as plt
import compbio.utils.memo as mem
def getct(n, seed = 0, forceRB = True,BW = False):
    random.seed(seed)
    ct = []
    for i in range(n):
        if forceRB:
            if i == 0:
                ct.append([1,0,0])
                continue
            elif i == 1:
                ct.append([0,0,1])
                continue
        ct.append([random.uniform(0,1),
                   random.uniform(0,1),
                   random.uniform(0,1)])
    return ct

def pmcolor(val,threshold = 0):
    if val < -1*threshold:
        return [1,0,0]
    elif val > threshold:
        return [0,0,1]
    else:
        return [0,0,0]

def blackbody(reset = False, **kwargs):
  '''Generate a colormap according to a logarithm of the blackbody spectrum. Colormap is transformed by an arctangent to place whites at .5 in a band with width determined by the [contrast] and then by taking the max with a gaussian peaked at 0,1 with width determined by [width]

kwargs: 
 reset
 contrast   [.64]   adjusts arctangent slope near zero.
 width      [.01]   adjust gaussian thresholding near endpoints.
 flip       [False] Hot areas are blue if flip is false.
 flip_ends  [False] Gaussian threshold swaps high and low colors.

For a usage example, see compbio/fun/ocean.py'''
  def setBB( **kwargs  ): 
    dstr = '2deg'
    import inspect
    import os
    import re

    contrast = kwargs.get('contrast', .64)
    width = kwargs.get('width', .01)
    flip_ends = kwargs.get('flip_ends', False)
    flip = kwargs.get('flip',False)
    thisdir= os.path.dirname(inspect.stack()[0][1])
    bbcols = [re.split(re.compile('\s+'),l) for l in \
              open(os.path.join(thisdir,'blackbody.tab')).readlines()
              if dstr in l]
    rgb = []
    for b in bbcols:
      rgb.append(array(b[7:10],float))
    
    rgb = array(rgb)
    npts = 512
    ntot = len(rgb)
    #DON'T TOUCH THE SCALE!!!!
    scl = 33.5
    
    #xvals = logspace(0.,scl,npts)/pow(10,scl) - .5
    xvals =  arctan(linspace(-contrast,contrast,npts))/pi*2 
    #scaling = (linspace(-1,1,npts)**3)*(linspace(-1,1,npts)**polypow)
    #width = .001
    scaling = exp( - ( 1-abs(linspace(-1,1,npts))) **2 / width)
    
    scaling *= array([-1 if x <0 else 1 for x in linspace(-1,1,npts)])
    if flip_ends:
      scaling *= -1

    d2 = [xvals[:], 
         scaling[:]]
    inds = argmax(np.abs(d2),0)
    xvals = array([d2[inds[i]][i] for i in range(len(inds))])
    
    xvals = (xvals * .5) + .5
    f0 = float(argmin(var(rgb,1))) / ntot
    k = -2.0 * (2. * f0  - 1) / (f0**2)  /2
    lspace = log(linspace(1, 1 + k,ntot)) / log(1 + k)
    #x0s = log(linspace(1,1 + scl,ntot))/log(1 + scl)
    
    vals =array([ interp(xvals,
                         lspace,
                         [e[i] for e in rgb])
                  for i in range(3)])
    midpoint = argmin(var(vals,0))
    
#    raise Exception()
    if flip: vals = vals[:,::-1]
    xs=linspace(0,1,npts)
    cdict = dict(
        red = [  (xs[i],  vals[0,i], vals[0,i]) for i in range(npts)],
        green = [ ( xs[i], vals[1,i], vals[1,i]) for i in range(npts)],
        blue = [  (xs[i], vals[2,i], vals[2,i]) for i in range(npts)])
    
    #out = zip(vals)
    cmap = matplotlib.colors.LinearSegmentedColormap('bb',cdict) 
    return cmap
  out = mem.getOrSet(setBB, reset = reset, **kwargs)
  return out
    
  
def imgcolor(arr ,normal = True, 
             #spectrum =True, 
             BW = False,
             alpha = False,
             color = [1,0,0]
             ):
    acopy = array(arr)
    if normal:
        acopy -= np.min(acopy)
        acopy /= np.max(acopy)
    

        
    if BW:
        cell = array([1,1,1])
        a3d = array(acopy[:,:,newaxis] * cell)
        return a3d
    elif alpha:
        cell = array([0,0,0,1])
        a4d = array(acopy[:,:,newaxis] * cell)
        cell2 = concatenate((color,[0]))
        a4d = a4d + cell2
        return a4d
    else:
        #Assume spectrum....
        a3d = array([[array([j,1,1]) for j in i] for i in acopy])
        rgb = mplcolors.hsv_to_rgb(a3d)
        return rgb
