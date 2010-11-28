import matplotlib.pyplot as plt
import numpy as np
from numpy import *
import matplotlib.patches as patches
  
def seismic(fig,
            yvecs_in, 
            colors = None,
            labels = None,
            scale_individual = False,
            xax = None,
            edgecolor = 'none',
            linewidth = 0,
            scale_none = 'False',
            subplot = None,
            axison = 'off',
            y_marked = None,
            continuous = True):
    
    
    #Seismic only shows absolute values... deal with it.
    yvecs = np.abs(yvecs_in)
    n = shape(yvecs)[0]
    nx = shape(yvecs)[1]
    yofs = arange(n)
   
    big_max = np.max(yvecs)
    all_maxes = np.max(yvecs,1)
    
    if not subplot:
        ax = fig.add_axes([0,0,1,1],frameon = False)
    else:
        ax = fig.add_subplot(subplot)
    if not xax:
        xax = arange(nx)

    xbounds = [min(xax),max(xax)]

    xax = array(xax)
    
    plt.axis(axison)

    #Just color everything blue if color is unspecified
    if not colors:
        colors = ['black' for i in range(n)]

    if not continuous:
        xinds = ((array((arange(2 * len(xax))),int)+1) /2)[:-1]
        yinds = ((array((arange(2 * len(xax))),int)/2))[:-1]
    else:
        xinds = array(arange(0,len(xax)),int)
        yinds = xinds

    nz_lambda = lambda x: x == 0 and 1 or x
    for i in range(n):
        yfrac = 1.0/2.5/n
        ypad = 1.1
        yofs =  (float(i)+.5) / n
        yvec = yvecs[i]

        

        if scale_individual:
            scl = 1.0/ nz_lambda(all_maxes[i])/ ypad * yfrac 
        elif not scale_none:
            scl = 1.0/nz_lambda(big_max)/ ypad *yfrac 
        else:
            scl = 1.0/ypad*yfrac
        yvec*=scl

        ax.fill_between(xax[xinds],yvec[yinds]+yofs,-1*yvec[yinds]+yofs,
                        edgecolor = edgecolor,
                        linewidth =linewidth,
                        facecolor = colors[i])


        if y_marked:
            


            xbounds = [ min(xax) - (max(xax) - min(xax))*.15,
                        max(xax) + (max(xax) - min(xax))*.15]
            xpts =[ min(xax) - (max(xax) - min(xax))*.1,
                    max(xax) + (max(xax) - min(xax))*.1]
            ypts = array([-y_marked,y_marked])* scl + yofs

            for x in xpts:
                p = patches.ConnectionPatch([x,ypts[0]],[x,ypts[1]],'data','data',
                                        arrowstyle = 'wedge,tail_width=1',
                                        edgecolor = 'black',
                                        facecolor = 'white',
                                        alpha = 1)
                ax.add_patch(p)
                ax.annotate(str(y_marked), [x, max(ypts)],xytext = [10,0],textcoords='offset points')

    ax.set_ylim([0,1])
    
    ax.set_xlim(xbounds)
    
    
