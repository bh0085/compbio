import netutils2 as nu2
import netutils as nu
from numpy import *
import numpy as np
import matplotlib.pyplot as plt
import utils.colors as mycolors

def show_binary(idx = 0):
    tsb = nu.expr_TS_binary(reset = 0)
    tsvals = nu.load_TS()

    net = nu2.get_net()
    tgs = net[1]
    tfs = net[0]

    f = plt.figure(0)
    f.clear()
    ax = f.add_subplot(111)


    
    for k in tsb.keys()[idx:]:
        
        my_tfs = tgs.get(k,[])
        ct = mycolors.getct(len(my_tfs))
        tgseries = tsvals[k]

        if not my_tfs: continue

        for i in range(len(my_tfs)):
            tf = my_tfs[i][0]
            
            series = tsvals.get(tf)
            if not series: continue
            binary = tsb.get(tf)
            #if not binary: 
            #    print 'no ts for ' + tg
            #    continue
            npts = len(binary)

            xax = tgseries
            cmap = equal(binary,0)[:,newaxis]*[1,0,0] + equal(binary,1)[:,newaxis]*[0,1,0]
        
            print my_tfs[i][1]
            ax.scatter(xax, series,  500, 
                       color = cmap,
                       alpha = my_tfs[i][1],
                       edgecolor = '0')
        break
    return
