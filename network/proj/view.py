import matplotlib.pyplot as plt
import netutils as nu
from numpy import *
import numpy as np
import compbio.utils.seismic as sm
import compbio.utils.plots as myplots

def view_in():
    na = nu.net_affinity()
    f = plt.figure(0)
    f.clear()
    ax = f.add_subplot(111)
    in_degree = sum(na,0)
    srt = argsort(in_degree)
    
    sm.seismic([in_degree[srt]],
               ax = ax)
    
    myplots.maketitle(ax,'In degree, sorted')
    
