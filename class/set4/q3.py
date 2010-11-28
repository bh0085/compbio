import re
import matplotlib.pyplot as plt
from numpy import *
import numpy as np
from utils import seismic, plots

def q3a(p = True):
    f = open('data/XPEHH.txt')
    f.readline()
    lines = f.readlines()
    print lines[-1]
    cols = ['ID', 'Chr','bp','cM','YvJC','CvJC','CvY']
    nc = len(cols)

    d = []
    for l in lines:
        m = re.search(re.compile('([^\s]+)\s*'*len(cols)),l)
        k = {}
        for i in range(nc):
            val = m.group(i+1)
            if i != 0: val = float(val)
            k[cols[i]] = val
        d.append(k)

    

    if p:

        labels = ['YvJC','CvJC','CvY']
        colors =  ['red','green','blue']
        f2 = plt.figure(2)
        f2.clear()
        seismic.seismic(f2, array(map(lambda x:[x['YvJC'],x['CvJC'],x['CvY']],d)).T, colors =colors, xax = map(lambda x:x['cM'],d), subplot = '111')
        
        plots.color_legend(f2, colors, labels)
    

    max_ehh = map(lambda x:np.max([x['YvJC'],x['CvJC'],x['CvY']]),d)
    mean_ehh = map(lambda x:np.mean([x['YvJC'],x['CvJC'],x['CvY']]),d)
    
    max_gt = nonzero(greater(max_ehh, 2))[0]
    mean_gt = nonzero(greater(mean_ehh,2))[0]
    
    print 'all pairwise:'
    print 'n bigger than ... 2'
    print 'max: ',len(max_gt)
    print 'mean: ',len(mean_gt)

    max_ehh = map(lambda x:np.max([x['CvJC'],x['CvY']]),d)
    mean_ehh = map(lambda x:np.mean([x['CvJC'],x['CvY']]),d)

    print 'Europe only'
    max_gt = nonzero(greater(max_ehh, 2))[0]
    mean_gt = nonzero(greater(mean_ehh,2))[0]
    

    print 'n bigger than ... 2'
    print 'max: ',len(max_gt)
    print 'mean: ',len(mean_gt)
    
    print 'max_gt:'
    print max_gt

    out = []
    for m in max_gt:
        out.append(d[m])

    return out


def q3b():
    f = open('data/Derived.txt')
    f.readline()
    lines = f.readlines()
    cols = ['ID', 'bp', 'cM', 'anc', 'der', 'C', 'Y', 'J']
    
    nc = len(cols)

    d = []
    for l in lines:
        m = re.search(re.compile('([^\s]+)\s*'*len(cols)),l)
        k = {}
        for i in range(nc):
            val = m.group(i+1)
            if not i in [0,3,4]: val = float(val)
            k[cols[i]] = val
        d.append(k)

    
    cfreq = array(map(lambda x: x['C'],d))/120
    yfreq = array(map(lambda x: x['Y'],d))/120
    jfreq = array(map(lambda x: x['J'],d))/180


    labels = ['C (Europe)','Y (Africa)','J (Asia)']
    colors =  ['red','green','blue']
    f2 = plt.figure(2)
    f2.clear()

    yvals = array([cfreq,yfreq,jfreq]) 
    ylim = seismic.seismic(f2,yvals , colors =colors, xax = map(lambda x:x['cM'],d),scale_none = True,axison='off', y_marked = .6)
    ax = f2.add_axes([0,0,1,1])
    plt.axis('off')
    ax.annotate('Derived allele frequencies on chromosome 5 (European)', [.05,.05],xytext=[.05,.98],textcoords='figure fraction',verticalalignment = 'top',family = 'serif', size = 'x-large')
    
    yvals2 = array(yvals)
    yvals2[nonzero(greater(yvals2,.6))] = .6

    #return

    plots.color_legend(f2, colors, labels)

    maxes = q3a(p = False)
    f2.clear()
    #ax = f2.add_subplot(111)

    idlam = lambda x:x['ID']
    long_ids = map(idlam,maxes) 
    these_ids = map(idlam,d)
    long_inds= [these_ids.index(i) for i in long_ids]
    dsub = [ d[i] for i in long_inds]

    
 

    seismic.seismic(f2,yvals[:,long_inds] , colors =colors,scale_none = True,subplot = '111', y_marked = .6, continuous = False)
    seismic.seismic(f2,yvals2[:,long_inds] , colors ='www',scale_none = True,subplot = '111',linewidth = 1,edgecolor = 'black',y_marked = .6,continuous = False)
    plots.color_legend(f2, colors, labels)
    plt.axis('off')
    ax = f2.add_axes([0,0,1,1])
    plt.axis('off')
    ax.annotate('Derived fraction for long haplotypes: Colors showed when derived fraction is .6',[.5,.5],xytext = [.05,.9],textcoords='figure fraction', color = 'black',family='serif',size='x-large')
    
    
    nbig = len(nonzero(greater(cfreq[long_inds],.6))[0])
    ax.annotate(str(nbig) + ' long haps for which frac derived > .6 in European pop',[.5,.5],xytext = [.95,.05],textcoords='figure fraction', color = 'black',family='serif',size='x-large',horizontalalignment='right')
