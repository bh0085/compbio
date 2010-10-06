import os 
import numpy as np
import itertools as it
import re   
import matplotlib.pyplot as plt
 

def show_m(net_number = 0):
    #A few functions defined here to be used later
    trgfun = lambda x: x[1]
    wtfun = lambda x:float( x[2] )
    tffun = lambda x: x[0]
    sigmafun = lambda x: 1 / (1 + np.exp(-x /1))

    prefix = 'predmodel/regressionwts/bRN'
    nwdata = open(os.path.join(prefix,'nw_'+str(net_number) + '.sif')).read()
    #Parse the list
    r = re.compile('^[ ]*(?P<target>[^\s]+)\t(?P<tf>[^\s]+)\t(?P<weight>[^\s]+)',re.M)
    matches = list(re.finditer(r,nwdata))    
    #Unsorted lists of tfs and targets
    targets =map(lambda x:x.group('target'),matches)
    tfs =    map(lambda x:x.group('tf'),matches)
    weights =map(lambda x:x.group('weight'),matches)
    
    #Concat the data for easier sorting
    cat = []
    for i in np.argsort(tfs):
        cat.append([tfs[i],targets[i],weights[i]])
    #Extract a dictionary with information for each target.
    trg_d = {}
    count = 0.0
    for k, g in it.groupby(sorted(cat,key = trgfun),key = trgfun):
        l = list(g)
        count += 1.0
        trg_d[k] = {'color': np.array([count, 0, 0]),
                    'tfs' : map(tffun,l),
                    'weights': map(wtfun,l)
                    }
    for k, v in trg_d.iteritems():
        v['color'] /= count
    trgnames =  trg_d.keys()
    
    import random
    trgcolors = {}
    for n in trgnames:
        c = [random.uniform(0,1),random.uniform(0,1),random.uniform(0,1)]
        trgcolors[n]= c

    #Extract a dictionary with information for each TF
    tf_d = {}
    for k, g in it.groupby(cat,key = lambda x: x[0]):
        l = list(g)
        tf_targets = map(lambda x: x[1],l)
        
        tf_d[k] = {'targets':tf_targets}

    tflens = map(lambda x: len(tf_d[x]['targets']), tf_d.keys())
    tfsort = np.argsort(tflens)
    

    #Display options:
    rtype = 'abs' #valid types, abs and val
    axis =  'polar' # valid types are rect and polar
    colortype = 'sign' # valid types are sign and identity

    #Configure the plotter
    fig = plt.figure(net_number,[8,8],facecolor = 'w')
    plt.clf()
    if axis == 'polar':
        ax = plt.axes(polar = True, frameon = False)
        ax.set_ylim([0,1.3])

    else:
        ax = plt.axes()
        ax.set_xlim([0,2 *np.pi])
        ax.set_ylim([0,1.5])
    inner_rad = .2
    plt.axis('off')


    #Choose subsets of interest for TF and Gene
    n = len(tfsort) #/5
    trg_subset = trgnames[::1]
    
    #Initialize the loop
    xs , ys, cs, rs,ts,ps = [],[],[],[],[],[]
    for i in range(n):
        j = tfsort[-1 -i]
        #mid = np.array(divmod(i, np.ceil(np.sqrt(n))) ) / np.ceil(np.sqrt(n))
        mid = np.array([np.cos(np.pi*2 * float(i) / n), 
                        np.sin(2*np.pi*float(i)/n)])
        rval = .5 + inner_rad
        xs.insert(0,mid[0]*rval)
        ys.insert(0,mid[1]*rval)
        ts.insert(0,2 * np.pi * i / n)
        ps.insert(0,rval)
        rs.insert(0,tflens[j])
        cs.insert(0,[0,0,1])
        tfname = tf_d.keys()[j]
        tf =  tf_d[tfname]
        these_w = []
        for trg in tf['targets']:

            #if len(trg_d[trg]['weights']) > 1:
            #    continue
            if not trg in trg_subset:
                continue
            
        
            trg_tf_idx = trg_d[trg]['tfs'].index(tfname)
            trg_tf_w = trg_d[trg]['weights'][trg_tf_idx]
            if colortype == 'identity':
                color = trgcolors[trg]
                #color = trg_d[trg]['color']
            elif colortype == 'sign':
                thr = .1
                if trg_tf_w > thr:
                    color = [0,1,0]
                elif trg_tf_w < -1* thr:
                    color = [1,0,0]
                else:
                    color = [0,0,0]

            tbump = 0.0
            if trg_tf_w > 0:
                tbump = .1
            these_w.append(trg_tf_w)
            if rtype == 'val':
                rval = sigmafun(trg_tf_w) + inner_rad
            else:
                rval = np.abs(trg_tf_w) + inner_rad
                
            xs.append(mid[0]*rval)
            ys.append(mid[1]*rval + tbump)
            ts.append(2 * np.pi * i / n +tbump)
            ps.append(rval)
            rs.append(10)
            cs.append(color)
            
        these_w = np.array(these_w)
        std = these_w.std()
        m = these_w.mean()
        lb = m - std
        ub = m - std
        

    #Plot
    ax.scatter(ts,ps,rs,cs,
               edgecolor = 'none')

    legs = (plt.plot(1,1,markeredgecolor = 'none',markerfacecolor = [1,0,0],
                     visible = 'false',linestyle = 'none',
                     marker = '.',markersize = 10,label = 'w-'),
            plt.plot(1,1,markeredgecolor ='none',markerfacecolor = [0,1,0],
                     visible = 'false',linestyle = 'none',
                     marker = '.',markersize = 10,label = 'w+'),
            plt.plot(1,1,markeredgecolor = 'none',markerfacecolor = [0,0,0],
                     visible = 'false',linestyle = 'none',
                     marker = '.',markersize = 10,label = 'w-'),
            plt.plot(1,1,markeredgecolor = 'none',markerfacecolor = [0,0,1],
                     visible = 'false',linestyle = 'none',
                     marker = '.',markersize = 30,label = 'TFH'),
            plt.plot(1,1,markeredgecolor = 'none',markerfacecolor = [0,0,1],
                     visible = 'false',linestyle = 'none',
                     marker = '.',markersize = 10,label = 'TFL'))

    ax.legend(legs,
              ('Enhancer','Inhibitor','No Effect','TF, High Out','TF Low Out'),
              numpoints = 1,markerscale = 2)    

    #raise Exception()
