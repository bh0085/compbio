import re
import matplotlib.pyplot as plt
from numpy import *
import numpy as np
from utils import seismic, plots, colors as mycolors
from numpy import logical_and as land

def q3a(p = True):
    f = open('data/XPEHH.txt')
    f.readline()
    lines = f.readlines()
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
        ax = seismic.seismic(f2, array(map(lambda x:[x['YvJC'],x['CvJC'],x['CvY']],d)).T, colors =colors, xax = map(lambda x:x['cM'],d),
                             xunits = 'cM',yunits = 'score',
                             subplot = '111',linewidth =5 )
        
        plots.color_legend(f2, colors, labels, ax = ax,pos = 0)
        plots.maketitle(ax, 'XP-EHH signal across chromosome 5 for pairwise comparisons','Plotted against cM')
    

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


def q3b(p = True):
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
    f2 = plt.figure(0)
    f2.clear()
    ax = f2.add_subplot(111)
    yvals = array([cfreq,yfreq,jfreq]) 
    
    ax = seismic.seismic(f2,[yvals[0,:]] , 
                           colors =colors, 
                           xax = map(lambda x:x['cM'],d),
                           scale_none = True,axison='off', 
                           y_marked = .6,
                           ax = ax,
                           xunits = 'cM',
                           yunits = 'freq')
    
    plots.maketitle(ax,'Derived allele frequency for Europe', '')

#plots.maketitle(ax,'Derived allele frequency on chromosome 5','hi')
    #plots.color_legend(f2, colors, labels, ax = ax)


    yvals2 = array(yvals)
    yvals2[nonzero(greater(yvals2,.6))] = .6


    maxes = q3a(p = False)
    f2.clear()
    ax = f2.add_subplot(111)

    idlam = lambda x:x['ID']
    long_ids = map(idlam,maxes) 
    these_ids = map(idlam,d)
    long_inds= [these_ids.index(i) for i in long_ids]
    dsub = [ d[i] for i in long_inds]

    
 
    
    seismic.seismic(f2,yvals[:,long_inds] , ax = ax, colors =colors,scale_none = True, y_marked = .6, continuous = False,yunits = 'freq',xunits = 'SNP index')
    seismic.seismic(f2,yvals2[:,long_inds] , colors ='www',scale_none = True,ax = ax,linewidth = 1,edgecolor = 'black',label_x = False, label_y= False,y_marked = .6,continuous = False)
    plots.color_legend(f2, colors, labels, ax = ax)
    plots.maketitle(ax,'Derived fraction for long haplotypes across all populations', 'Colors showed when derived fraction is >.6')
    
    nbig = len(nonzero(greater(cfreq[long_inds],.6))[0])
    ax.annotate(str(nbig) + ' long haps for which frac derived > .6 in European pop',[.5,.5],xytext = [.95,.05],textcoords='figure fraction', color = 'black',family='serif',size='x-large',horizontalalignment='right')
    
    long_inds = array(long_inds)
    idxs_passing = array(long_inds[greater(cfreq[long_inds],.6)],int)
    return idxs_passing
    

def q3c():
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

    dout = list(d)

    
    cfreq = array(map(lambda x: x['C'],d))/120
    yfreq = array(map(lambda x: x['Y'],d))/120
    jfreq = array(map(lambda x: x['J'],d))/180

    pair_y = array([cfreq,yfreq])
    pair_j = array([cfreq,jfreq])
    
    f = plt.figure(0)

    fs = zeros((3,len(cfreq)))
    for i in range(2):
        pairing = [pair_y,pair_j][i]
        p = mean(pairing,0)
        d = cfreq - p

        fs[i,:]  = power(d,2)/(p*(1-p))
        
    fs[isnan(fs)] = 0 

    fs[2,:] = np.mean(fs[0:2,:],0)
    f.clear()
    ct = mycolors.getct(3)

    ax = seismic.seismic(f, fs, subplot = 111,
                     colors = ct,continuous = False,
                    yunits = 'Fst'
                    ,xunits = 'position')
        
    plots.color_legend(f,ct,[r'European-Asian $F_{ST}$',
                      r'European-African $F_{ST}$',
                      r'Average $F_{ST}$.'])
    plots.maketitle(ax,'$F_{ST}$ Across chromosome 5')
    

    idxs_passing = q3b()
    f = plt.figure(0)
    f.clear()
    #ax2 = f.add_subplot(111,frameon = False)

    ax2 = seismic.seismic(f, fs[:,idxs_passing],subplot = 111,
                     colors = ct,continuous = False, scale_none = True,
                    yunits = 'Fst',
                    xunits = 'position',axison = 'off')
        
    fs2 = array(fs)

    previous_tests = zeros(len(fs2[0]))
    previous_tests[idxs_passing] = 1
    print 'NBIG:   ', len(nonzero(land(greater(fs2[2,:],.6), previous_tests))[0])

    fs2[greater(fs2,.6)] = .6
    
    ax2 = seismic.seismic(f, fs2[:,idxs_passing], subplot = 111, scale_none = True,
                     colors ='www',continuous = False,label_x = False, label_y = False,
                    yunits = 'Fst' 
                    ,xunits = 'position')

    plots.color_legend(f,ct,[r'European-Asian $F_{ST}$',
                      r'European-African $F_{ST}$',
                      r'Average $F_{ST}$.'], ax = ax2)
    plots.maketitle(ax2,'$F_{ST}$ At spots having long haplotype and derived allele frequency signals', 'Colors indicate values greater than .6.\n' + str(len(nonzero(land(equal(fs2[2,:],.6),previous_tests))[0])) + ' sites for which all thresholds are passed')

    return dout, nonzero(land(equal(fs2[2,:],.6),previous_tests))[0]

def q3d(d,final_idxs):
    cons = open('data/phastConsElements.txt')
    cons.readline()
    cons = cons.readlines()
    genes = open('data/genes.gff')
    genes.readline()
    genes = genes.readlines()
    

    gcols = ['chr','source','category','start','end','score','strand','phase','name']
    conscols = ['chr','start','stop','LOD']
    

    nc = len(conscols)

    consd = []
    for l in cons:
        m = re.search(re.compile('([^\s]+)\s*'*len(conscols)),l)
        k = {}
        for i in range(nc):
            val = m.group(i+1)
            if not i in [3]: val = float(val)
            k[conscols[i]] = val
        consd.append(k)

    nc = len(gcols)
    genesd = []
    for l in genes:
        m = re.search(re.compile('([^\s]+)\s*'*len(gcols)),l)
        k = {}
        for i in range(nc):
            val = m.group(i+1)
            if not i in [1,2,5,6,7,8]: val = float(val)
            k[gcols[i]] = val
        genesd.append(k)

    g5 = []
    for g in genesd:
        if g['chr'] ==5:
            g5.append(g)

            
    excount = 0
    intcount = 0
    conscount = 0 
    for snp in d:
        bp = snp['bp']
        cmatches = []
        gmatches = []
        for g in g5:
            if g['category'] == 'exon' or g['category'] == 'intron':
                if bp > g['start'] and bp < g['end']:
                    gmatches.append(g)
        for c in consd:
            if bp > c['start'] and bp < c['stop']:
                cmatches.append(c)        

        if cmatches:
            conscount += 1
        if gmatches:
            for g in gmatches:
                if g['category'] == 'intron':
                    intcount +=1 
                    break
            for g in gmatches:
                if g['category'] == 'exon':
                    excount +=1
                    break

    print len(d)
    print 'excount: ',excount
    print 'intcount: ',intcount
    print 'conscount: ',conscount

    for idx in final_idxs:
        bp = d[idx]['bp']
        
        gmatches = []
        cmatches = []
        for g in g5:
            if g['category'] == 'exon' or g['category'] == 'intron':
                if bp > g['start'] and bp < g['end']:
                    gmatches.append(g)
        for c in consd:
            if bp > c['start'] and bp < c['stop']:
                cmatches.append(c)
        import pprint
        pprint.pprint( d[idx] )
        pprint.pprint( gmatches )
        print

        for g in gmatches:
            if g['category'] == 'exon': print "EXON"
            if g['category'] == 'intron': print "INTRON"

        if cmatches:
            pprint.pprint(cmatches)
            for g in gmatches:
                print
                print g['start'], bp, g['end'],g['name']
        else:
            print('no conservation matches for SNP')
