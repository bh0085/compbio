import matplotlib.pyplot as plt
import numpy as np
from numpy import *
import matplotlib.patches as patches
  


def seismic(yvecs_in, 
            colors = None,
            labels = None,
            scale_individual = False,
            xax = None,
            edgecolor = 'none',
            linewidth = 0,
            scale_none = False,
            subplot = None,
            axison = 'off',
            y_marked = None,
            mark_scale = True,
            continuous = True,
            ax = None,
            yunits = '',
            xunits = '',
            label_y = True,
            label_x = True,
            labelspace= False,
            labelinds =None,
            fig = None,
            stacked = False,
            ymarkpts = None,
            xmarkpts = None,
            xmark_colors = None,
            right_label_on = False):

        

    #Seismic only shows absolute values... deal with it.
    yvecs = np.abs(yvecs_in)
    n = shape(yvecs)[0]
    nx = shape(yvecs)[1]
    yofs = arange(n)
   
    big_max = np.max(yvecs)
    all_maxes = np.max(yvecs,1)
    
    if labelinds == None:
        labelinds = arange(n)

    if not ax:
        if not fig: 
            fig = plt.figure(0)
            fig.clear()

        if not subplot:
            ax = fig.add_axes([0,0,1,1],frameon = False)
        else:
            ax = fig.add_subplot(subplot)
    if xax == None:
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
    ysofar = zeros(nx)
    max_y = np.max(yvecs)
    max_sum_y = np.max(np.sum(yvecs,0))
    for i in range(n):
        if not labelspace: fracmax = 1
        else: fracmax = .75

        yfrac = 1.0/2.0/(n) * fracmax
        ypad = 1.1
        yvec = yvecs[i]

        if not stacked:
            yofs =  (float(i) + .5) /(n) *fracmax

            if scale_individual:
                scl = 1.0/ nz_lambda(all_maxes[i])/ ypad * yfrac 
            elif not scale_none:
                scl = 1.0/nz_lambda(big_max)/ ypad *yfrac 
            else:
                scl = 1.0/ypad*yfrac

            yup = yvec[yinds]*scl+yofs
            ydown = -1*yvec[yinds]*scl+yofs
            ax.fill_between(xax[xinds],yup,ydown,
                            edgecolor = edgecolor,
                            linewidth =linewidth,
                            facecolor = colors[i])
            if xmarkpts != None:
              if not shape(xmarkpts): these_marks = [xmarkpts]
              elif not shape(xmarkpts[i]): these_marks = [xmarkpts[i]]
              else: these_marks = xmarkpts[i]
              if not xmark_colors:
                xmark_colors = ['black'] * len(these_marks)
              for idx,xm in enumerate(these_marks):
                xmark_ys = [interp( xmarkpts[i],xax[xinds], yup),
                            interp(xmarkpts[i], xax[xinds], ydown)]

                ax.plot([xm] *2 , xmark_ys, color =xmark_colors[idx], linewidth = 2)
                print xmark_ys, [xm] * 2
                print xax[xinds]
              
        else:
            if i == 0: ydown = array(ysofar)
            yofs = 0
            scl = 1.0/max_sum_y
            yup = (yvec + ysofar)*scl +yofs
            ybelow = ysofar *scl + yofs
            
            
            ax.fill_between(xax[xinds],yup[yinds],ybelow[yinds],
                            edgecolor = edgecolor,
                            linewidth =linewidth,
                            facecolor = colors[i]
                            )
            ysofar += yvec
            ymarkpts = [min(ysofar),max(ysofar)]

        if not stacked:
            if i in labelinds: dolabel = True
            else: dolabel = False
        elif i == n -1:    dolabel = True
        else:  dolabel = False

        if dolabel:

            if y_marked:
                this_ymark = y_marked
            else:
                this_ymark = max(yvec[yinds])
            
            xbounds = [ min(xax) - (max(xax) - min(xax))*.15,
                        max(xax) + (max(xax) - min(xax))*.15]
            xpts =[ min(xax) - (max(xax) - min(xax))*.1,
                    max(xax) + (max(xax) - min(xax))*.1]
            if ymarkpts == None:
                ypts = array([-this_ymark,this_ymark])* scl + yofs
            else:
                ypts = array(ymarkpts)*scl + yofs
            
            dx = max(xax) - min(xax)
            dy = max(yup) - min(ydown)
            
            if label_y:
                for idx, elt  in enumerate(zip(xpts, ('right','left'))):
                    x, ha = elt
                    if idx == 1 and right_label_on == False:
                        continue
                    p = patches.ConnectionPatch([x,ypts[0]],[x,ypts[1]],'data','data',
                                            arrowstyle = '<->',
                                            edgecolor = 'black',
                                            facecolor = 'white',
                                            alpha = .5)
                    ax.add_patch(p)
                    ystr ='{0:2.2g}'.format(this_ymark)
                    if yunits != '':
                        ystr+=str('\n(')+str(yunits)+')'
                    ax.annotate(ystr, 
                                [x, max(ypts)],
                                xytext = [10* (-1 if ha == 'right' else 1),0],
                                textcoords='offset points',
                                verticalalignment = 'top',
                                ha = ha)
            
            if i == 0 and label_x:
            
                p = patches.ConnectionPatch([max(xax),min(ydown) - dy * .2],
                                             [min(xax),min(ydown) - dy * .2],
                                             'data','data',
                                            arrowstyle = '<->',
                                            edgecolor = 'black',
                                            facecolor = 'white',
                                            alpha = .5)
                ax.add_patch(p)
            
                ax.annotate(xunits,[min(xax) + .5*dx, min(ydown) - dy*.3],xycoords = 'data',
                            xytext = [.5,0],textcoords ='offset pixels', 
                            horizontalalignment = 'center',verticalalignment = 'top')
            
                ax.annotate(str(min(xax)),[min(xax), min(ydown) - dy*.3],xycoords = 'data',
                            xytext = [10,-5],textcoords = 'offset pixels', 
                            verticalalignment = 'top',horizontalalignment = 'center')
            
                ax.annotate(str(max(xax)),[max(xax), min(ydown) - dy*.3],xycoords = 'data',
                            xytext = [-10,-5],textcoords = 'offset pixels', 
                            verticalalignment = 'top', horizontalalignment = 'center')            
            
            if labels: 
                this_label = labels[i]
                ax.annotate(this_label, [min(xpts),mean(ypts)], 
                            size = 'x-large', textcoords = 'offset pixels',
                            xytext = [-10,0],color = colors[i],
                            horizontalalignment = 'right',
                            verticalalignment = 'center')

    ax.set_ylim([0,1])
    
    ax.set_xlim(xbounds)
    return ax
    
    
