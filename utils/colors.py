import random
from numpy import *
import numpy as np
import matplotlib.colors as mplcolors
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
