import random
from numpy import *
import numpy as np
import matplotlib.colors as mplcolors
import matplotlib
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

def blackbody(contrast = 1.5, polypow = 5):
  #def setBB():
  #  bbpath = 
  dstr = '2deg'
  import inspect
  import os
  import re
  thisdir= os.path.dirname(inspect.stack()[0][1])
  bbcols = [re.split(re.compile('\s+'),l) for l in \
            open(os.path.join(thisdir,'blackbody.tab')).readlines()
            if dstr in l]
  rgb = []
  for b in bbcols:
    rgb.append(array(b[7:10],float))
  
  rgb = array(rgb)
  npts = 256
  ntot = len(rgb)
  #DON'T TOUCH THE SCALE!!!!
  scl = 34.5
  
  #xvals = logspace(0.,scl,npts)/pow(10,scl) - .5
  xvals =  arctan(linspace(-contrast,contrast,npts))/pi*2 
  scaling = (linspace(-1,1,npts)**3)*(linspace(-1,1,npts)**polypow)
  d2 = [xvals[:], 
       scaling[:]]
  inds = argmax(np.abs(d2),0)
  xvals = array([d2[inds[i]][i] for i in range(len(inds))])
  
  
  xvals = (xvals * .5) + .5

  matplotlib.pyplot.plot(xvals)
  matplotlib.pyplot.plot(xvals[::-1])
#.5
  #xvals =  xvals - min(xvals)
  #xvals /= max(xvals)

  #print xvals[499]
  #print xvals[500]

  #print xvals
  x0s = log(linspace(1,scl,ntot))/log(scl)

  vals =array([ interp(xvals,
                       ,x0s,
                       [e[i] for e in rgb])
          for i in range(3)])

  return vals
  xs=linspace(0,1,npts)

  cdict = dict(
      red = [  (xs[i],  vals[0,i], vals[0,i]) for i in range(npts)],
      green = [ ( xs[i], vals[1,i], vals[1,i]) for i in range(npts)],
      blue = [  (xs[i], vals[2,i], vals[2,i]) for i in range(npts)])
  
  #out = zip(vals)
  cmap = matplotlib.colors.LinearSegmentedColormap('bb',cdict) 
  
  return cmap
    
  
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
