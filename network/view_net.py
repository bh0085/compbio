import netutils as nu
import netsvd as ns
import matplotlib.pyplot as plt
from numpy import *
import numpy as np

def view():
    arr = reshape(arange(100.),(10,10))
    view_array(arr)

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
            r = abs(e) * rscl
            c = [int(e > 0), 0 , int(e < 0)]
            xs.append(x)
            ys.append(y)
            rs.append(r)
            cs.append(c)

    fig = plt.figure(1)
    fig.clear()
    sp = fig.add_subplot(221,frame_on = False); sp.set_axis_off()
    print xs, ys,rs,cs
    sp.scatter(xs,ys,rs,cs)
    
