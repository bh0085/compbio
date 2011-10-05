import cb.utils.plots as myplots
from numpy import *
import numpy as np
import os

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def load(res = 25):
    if res == 25: fpath = '/data/brain_atlas/AtlasAnnotation25.sva'
    else: raise Exception()

    print 'path: ', fpath
    size =  os.path.getsize(fpath)
    n = 10000
    skip = size / n 
    
    f = open(fpath)
    f.readline()
    f.readline()
    coords = []
    evals = []
    while( len(coords) < n):
        
        f.seek(skip, 1)
        l0 = f.readline()
        l = f.readline()
        if( l == ''): break;
        lvals =[float(v) for v in l.split(',')]
        coords.append(tuple(lvals[0:3]))
        evals.append(lvals[3])
        

    print 'len: ', (len(coords))
    fig = myplots.fignum(1,(8,6))
    
    ax = fig.add_subplot(111, projection='3d')

    xyvals = array([[x[0],x[1], x[2]] for x in coords])
    evals = array(evals)
    evals = (evals / np.max(evals) )[:,newaxis] * array([1.,0,0])
    print shape(xyvals)
    ax.scatter(xyvals[:,0], xyvals[:,1], xyvals[:,2], 
               s = 5,
               edgecolor = 'none',
               facecolor = evals)

    #raise Exception()

    path = myplots.figpath('brainmap_spatial_coords_{0}'.format(res))
    fig.savefig(path)
        


    return coords
