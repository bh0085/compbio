import bsort as bs
from numpy import *
import cb.utils.plots as myplots
import numpy as np
import matplotlib.pyplot as plt

dfile = '/data/misc/CorrelationPlot.txt'
data = open(dfile).readlines()
arr = array([[float(e) for e in l.split('\t')[1:]] 
             for l in data[1:]])

rows =reshape(array([[e for e in l.split('\t')[:1]] 
              for l in data[1:]]),-1)
cols =reshape(array([[e for e in l.split('\t')[:]] 
              for l in data[:1]]),-1)
#arr = arr * greater(arr, 0)
arr = arr/np.max(arr,0)


def run(meth = 'moment'):
    out,srts = bs.run0(arr = arr, itr = 2, meth = meth)
    f = myplots.fignum(3,(12,6))
    ax = f.add_subplot(111)

    csrts = [s for s in srts if len(s) == len(cols)][0]
    rsrts = [s for s in srts if len(s) == len(rows)][0]
    cprint = [rows[rs] for rs in rsrts]
    rprint = [cols[cs] for cs in csrts]


    im = ax.imshow(out,
              interpolation= 'nearest',
              cmap = plt.get_cmap('OrRd'),
              )

        #flip the rows and columns... looks better.   
    ax.set_xticks(arange(len(cols))+.25)
    ax.set_yticks(arange(len(rows))+.25)

    ax.set_yticklabels([e for  e in cprint])
    ax.set_xticklabels(rprint)

    print 'rows: \n{0}'.format(', '.join([e.strip() for e in rprint]))
    print
    print 'cols: \n{0}'.format(', '.join([e.strip() for e in cprint]))

    plt.colorbar(im)
    
    f.savefig(myplots.figpath('correlation_plot_2_4_{0}.pdf')
              .format(meth))
    return
