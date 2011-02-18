#A sample QT GUI within which to view matrices.

import netutils as nu
import netsvd as ns
import matplotlib.pyplot as plt
from numpy import *

import sys, os, random

import compbio.utils.colors as colors
import matplotlib
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure



def view_array(M):
    
    m = M.shape[0]
    n = M.shape[1]
    

    xs,ys,rs,cs = [], [], [], []
    max_r = 50^2
    rscl = max_r / (abs(M).max())
    
    
    for i in range(m):
        for j in range(n):
            y = i
            x = j
            e = M[i,j]
            r = 100
            c =colors.pmcolor(e,3)
            xs.append(x)
            ys.append(y)
            rs.append(r)
            cs.append(c)

    fig = plt.figure(1)
    fig.clear()
    sp = fig.add_subplot(221,frame_on = False); sp.set_axis_off()
    sp.scatter(xs,ys,rs,cs)
    plt.show()

