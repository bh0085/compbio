import os
import numpy as np
import re
import itertools as it
import scipy as sp
import scipy.linalg as lin
import matplotlib.pyplot as plt
from scipy import sparse
from numpy import *
import pickle
import subprocess
import random
import mlpy
import compbio.utils.colors as mycolors
import netwriter as nw
import inspect
import pickle
import inspect

default_name = 'unsup'

def parse_net(name = default_name,number = 0, reset = 0):
    
    net_dir = os.path.abspath(os.path.dirname(inspect.getfile(inspect.currentframe())))
    if not reset:
        net,sxs = nw.rn2(name = name)
        if not sxs: raise Exception('error reading ' + name)
    else:
        if name == 'unsup':
            fpath = os.path.join(net_dir,'patrick/unsup_patrick.txt')
        else:
            fpath = os.path.join(net_dir,'patrick/logistic_0.6.txt')


        TS = load_TS( reset = mod(reset,2))
        CL = load_CL( reset = mod(reset,2))
        nwdata = open(fpath).read()
    #A few functions defined here to be used later
        trgfun = lambda x: x[1]
        wtfun = lambda x:float( x[2] )
        tffun = lambda x: x[0]
        sigmafun = lambda x: 1 / (1 + np.exp(-x /1))

        r = re.compile('^[ ]*(?P<tf>\S+)\s+(?P<target>\S+)\s+(?P<weight>\S+)'
                   ,re.M)
        matches = list(re.finditer(r,nwdata))    
    #Unsorted lists of tfs and targets
        targets =map(lambda x:x.group('target'),matches)
        tfs =    map(lambda x:x.group('tf'),matches)
        weights =map(lambda x:x.group('weight'),matches)
    
    #Concat the data for easier sorting
        cat = []
        for i in np.argsort(tfs):
            if TS.has_key(tfs[i]) and CL.has_key(targets[i]):
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


    #Extract a dictionary with information for each TF
        tf_d = {}
        for k, g in it.groupby(cat,key = lambda x: x[0]):
            l = list(g)
            tf_targets = map(lambda x: x[1],l)
        
            tf_d[k] = {'targets':map(trgfun,l),
                       'weights':map(wtfun,l)}

        net =  (trg_d, tf_d)
        nw.wn2( name,net, hardcopy = True)
    return net


def parse_TS():
    f = open('TC.geneexp').read()
    elts =f.split('\n')
    seqdict = {}
    for e in elts:
        matches = list(re.finditer(re.compile('([^\s]+)'), e))
        if not len(matches): continue
        name = matches[0].group(1)
        seqdict[name] = []
        for i in matches[1:]:
            seqdict[name].append(float(i.group(1)))

    pickle.dump(seqdict,open('TC.pickle','w'))

def parse_CL():
    f = open('CL.geneexp').read()
    elts =f.split('\n')
    seqdict = {}
    for e in elts:
        matches = list(re.finditer(re.compile('([^\s]+)'), e))
        if not len(matches): continue
        name = matches[0].group(1)
        seqdict[name] = []
        for i in matches[1:]:
            seqdict[name].append(float(i.group(1)))

    pickle.dump(seqdict,open('CL.pickle','w'))

def load_CL(reset = 0):
    hardcopy = True
    net_dir = os.path.abspath(os.path.dirname(inspect.getfile(inspect.currentframe())))

    if not reset:
        #no reason to use name... only one cl is available
        out, sxs = nw.rn2( default_name, hardcopy = hardcopy, np = False)
        if not sxs: raise Exception()
    else:
        out = pickle.load(open(os.path.join(net_dir, 'CL.pickle')))
        nw.wn2( default_name, out,  hardcopy = hardcopy, np = False)

    return out

def load_TS(reset = 0):
    hardcopy = True
    net_dir = os.path.abspath(os.path.dirname(inspect.getfile(inspect.currentframe())))
    if not reset:
        #no reason to use name... only one cl is available
        out, sxs = nw.rn2(default_name, hardcopy = hardcopy, np = False)
        if not sxs: raise Exception()
    else:
        out = pickle.load(open(os.path.join(net_dir,'TC.pickle')))
        nw.wn2(default_name , out, hardcopy = hardcopy, np = False)

    return out




#RESETS:
#reset value specifies both whether a function
#will recompute its value (reset > 0)
#and whether it will ask its children to recompute (reset mod 2)
def net_square_affinity(name = default_name, reset = 0):

    hardcopy = True
    if not reset:
        net,sxs= nw.rn2(name,hardcopy = hardcopy, np = np)
        if not sxs: raise Exception()
    else:
        nw.claim_reset()
        trgs, tfs = parse_net(name)
        k0 = trgs.keys()
        k1 = tfs.keys()
    

        all_keys = k0
        all_keys.extend(k1)
        all_keys = sort(all_keys)
        uk = []
        for g,v in it.groupby(all_keys): uk.append(g)
        nk = size(uk)
        N =zeros((nk,nk),float32)

        tf_imap = {}
        for k in tfs.keys():
            tf_imap[k] = uk.index(k)

        for k,t in trgs.iteritems():
            i = uk.index(k)
            for ktf in t['tfs']:
                j = tf_imap[ktf]
                N[i,j] = 1
        

        kidxs = {}
        for i in range(len(uk)):
            kidxs[uk[i]] = i
            
        net_sq_keyidxs(name = name, write = kidxs)


        net = N
        nw.wn2(name, net,hardcopy = hardcopy,np = True)
    return net

def net_affinity(name = default_name ,reset = 0):
    hardcopy = True
    if not reset:
        net, sxs = nw.rn2(name, hardcopy = hardcopy, np = True)
        if not sxs: raise Exception()
    else:
        nw.claim_reset()
        trgs, tfs = parse_net(name, reset = mod(reset,2))
        k0 = trgs.keys()
        k1 = tfs.keys()
    
        kmap0 = {}
        kmap1 = {}
        for i in range(len(k0)):
            kmap0[k0[i]] = i
        for i in range(len(k1)):
            kmap1[k1[i]] = i

    #Now, for a first shot, lets construct an n*m matrix
    #with possible interactions for each gene appearing in 
    #the network        
    #prune the network for testing purposes.
        
        net_trg_keyidxs(name = name, write = kmap0)
        net_tf_keyidxs(name = name, write = kmap1)

        m = len(kmap0.keys())
        
        n = len(kmap1.keys())
        N = np.zeros([m,n],float32)
        for k_trg,v_trg in trgs.items():
            i = kmap0[k_trg]
            for idx in range(len(v_trg['tfs'])):
                k_tf = v_trg['tfs'][idx]
                w_tf = v_trg['weights'][idx]
                j = kmap1[k_tf]
                N[i][j] = w_tf
                print w_tf

        net = N
        nw.writenet(name, net,hardcopy = hardcopy, np = True)
    return net

def net_jcf(name = default_name, reset = 0):
    hardcopy = False
    try:
        if reset: raise Exception('compute')
        return nw.readnet(name,hardcopy = hardcopy)
    except Exception as e:
        if e.args[0] != 'compute': raise Exception()
        nw.claim_reset()

        square_aff = net_square_affinity(name, reset = mod(reset,2))
        import pymat
        command = 'pymatjordan'
        arr = square_aff
        
        inp = {'array':arr}
        output = pymat.runmat(command, inp)
        nw.writenet(name, output, hardcopy = hardcopy)
        return output

    




def net_genegene(name = default_name, reset = 0):
    hardcopy = True

    if not reset:
        out, sxs = nw.rn2(name,hardcopy = hardcopy, np = True)
        if not sxs: raise Exception()
    else:
        nw.claim_reset()
        sq = net_square_affinity(name, reset = mod(reset,2))
        genegene = dot(sq,sq.T)
        out = genegene
        nw.wn2(name,genegene, hardcopy = hardcopy, np =True)
    return out



#Write keys of targets to a file when an affinity matrix is created
def net_sq_keyidxs(name = default_name, write = None ):
    if write != None:
        nw.wn2(name, write, hardcopy = True, np = False)
        idxs = write
    else:
        idxs,sxs = nw.rn2(name, hardcopy = True, np = False)
        if not sxs: raise Exception()
    return idxs


#Write keys of targets to a file when an affinity matrix is created
def net_trg_keyidxs(name = default_name, write = None ):
    if write != None:
        nw.wn2(name, write, hardcopy = True, np = False)
        idxs = write
    else:
        idxs,sxs = nw.rn2(name, hardcopy = True, np = False)
        if not sxs: raise Exception()
    return idxs

#Write keys of tfs to a file when an affinity matrix is created
def net_tf_keyidxs(name =default_name, write = None):            
    if write != None:
        nw.wn2(name, write, hardcopy = True, np = False)
        idxs = write
    else:
        idxs,sxs = nw.rn2(name, hardcopy = True, np = False)
        if not sxs: raise Exception()
    return idxs


def net_genegene_norm(name = default_name, reset = 0):
    hardcopy = True
    if not reset:
        out, sxs =  nw.rn2(name,hardcopy = hardcopy,np = True)
        if not sxs: raise Exception()
    else:
        nw.claim_reset()
        sq0 =array(net_square_affinity(name, reset = mod(reset,2))) 
        sq1 =array( sq0.T)
        n = shape(sq0)[0]

        for i in range(n):
            sq0[:,i] /= max(1,np.sum(sq0[:,i]))
            sq1[:,i] /= max(1,sp.sum(sq1[:,i]))
        genegene = dot(sq0,sq1)
        nw.writenet(name,genegene, hardcopy = hardcopy,np = True)
        out = genegene
    return out



def net_hcluster_genegene_norm(name = default_name, reset = 0, n = 200,
                              cosine = True,
                              hier = True,
                              cut = .5):
    hardcopy = True
    try:
        if reset: raise Exception('compute')
        return nw.readnet(name,hardcopy =hardcopy)
    except Exception as e:
        if e.args[0] != 'compute': raise Exception()
        nw.claim_reset()
        sqa = net_genegene_norm(name, reset = mod(reset,2))
        gg = sqa
        sub = matrix(gg[0:n,0:n])


        if cosine:
            for view in sub:
                view = view / sum(view)
        
        h = mlpy.HCluster(link = 'complete')
        comp = h.compute(sub)
        clustered = h.cut(cut)
   
        value = {'HCluster':h,'clusters':clustered}
        nw.writenet(name,value, hardcopy = hardcopy)
        return value
        
def net_kcluster_genegene_norm(name = default_name, reset = 0, n = 200, 
                               k = 5,
                               cosine = True):
    hardcopy = True
    try:
        if reset: raise Exception('compute')
        return nw.readnet(name,hardcopy =hardcopy)
    except Exception as e:
        if e.args[0] != 'compute': raise Exception()
        nw.claim_reset()
        sqa = net_genegene_norm(name, reset = mod(reset,2))
        gg = sqa
        sub = matrix(gg[0:n,0:n])


        if cosine:
            for view in sub:
                view = view / sum(view)


        
        kmeans = mlpy.Kmeans(k)
        clustered = kmeans.compute(sub)
        means = kmeans.means

        value = {'KMeans':kmeans,'clusters':clustered}
        nw.writenet(name,value, hardcopy = hardcopy)
        return value
        




def net_blanket( name = default_name , reset = 0, order = 1):


    saff = net_square_affinity(name,reset = mod(reset,2))

    m = saff
    n = shape(m)[0]
    for i in range(n):
        m[i,i]  = 1
        

    m2 = m * m.T
    lower100 = m2[0:500,0:500].todense()
    f = plt.figure(0)
    f.clear()
    ax = f.add_axes([0,0,1,1], aspect = 'auto',frameon=False)
    ax.imshow(lower100,interpolation = 'nearest')

    
    return


def expr_TS_binary(reset = 0 , 
                   name = default_name,
                   hardcopy = True ):
    if not reset:
        out, sxs =  nw.rn2(name,hardcopy = hardcopy)
        if not sxs: raise Exception()
    else:
        nw.claim_reset()        
        ts = load_TS()
        ts_mixtures = {}
        count = 0
        
        nc = len(ts[ts.keys()[0]])
        all_states = {}
        mixes  = {}
        
        for k in ts.keys():
            if mod(count,100) == 0 : print count
            count += 1
            t = ts[k]
            #ts_mixtures contains the probabilities of a classification
            #at each timepoint (0 = off, 1 = onn)
            mixture = expr_gmm_onoff(t)
            cutoff = .25 #demand 25% certainty for a label
            mix0 = squeeze(mixture[:,0])
            mix1 = squeeze(mixture[:,1])
            
            states = zeros(nc) 
            states[nonzero(greater(mix1, mix0 + cutoff))[0]] = 1
            states[nonzero(greater(mix0, mix1 + cutoff))[0]] = -1
            mix = (mix0,mix1)
            mixes[k] = mix
            all_states[k] = states
            
            if mod(count, 10) == 0: print count

        nw.wn2( name, (all_states,mixes), hardcopy = hardcopy)
        out= (all_states,mixes)
    return out
def expr_CL_binary(reset = 0 , 
                   name = default_name,
                   hardcopy = True ):
    try:
        if reset: raise Exception('compute')
        return nw.readnet(name,hardcopy = hardcopy)
    except Exception as e:
        if e.args[0] != 'compute': raise Exception()
        nw.nw.claim_reset()        
        ts = load_CL()
        ts_mixtures = {}
        count = 0
        
        nc = len(ts[ts.keys()[0]])
        all_states = {}
        
        for k in ts.keys():
            count += 1
            t = ts[k]
            #ts_mixtures contains the probabilities of a classification
            #at each timepoint (0 = off, 1 = onn)
            mixture = expr_gmm_onoff(t)
            cutoff = .25 #demand 25% certainty for a label
            mix0 = squeeze(mixture[:,0])
            mix1 = squeeze(mixture[:,1])
            
            states = zeros(nc) + 2
            states[nonzero(greater(mix1, mix0 + cutoff))[0]] = 1
            states[nonzero(greater(mix0, mix1 + cutoff))[0]] = 0
            all_states[k] = states

        nw.writenet( name, all_states, hardcopy = hardcopy)
        return all_states



def expr_view():
    ts = load_TS()
    
    f = plt.figure(1)
    f.clear()
    ax = f.add_axes([0,0,1,1])
    for k in tfkeys:
        
        expr_gmm_onoff(ts[k],draw = True)
        expr_gmm_onoff(ts[k],log_expr = False,draw =True, fig = 2)

        inp = raw_input('next')
        if inp == 'd':
            break
        
def expr_gmm_onoff(expr_in,
                   log_expr = False, 
                   fig = 1,
                   draw = False):
    
    expr = (array(expr_in))
    dev = std(expr)
    if log_expr:
        expr = log(expr + dev)
    n = len(expr)
    expr_array = zeros((n,1))
    for i in range(n):
        expr_array[i] = expr[i]
    expr = expr_array
        
    from scikits.learn import gmm
    
    #demand seperation of max from alternate hypotheses by e/2
    cmin_diff = log(e*1.5)
    #cmin_diff = .0001
    k = 2

    G = gmm.GMM(n_states = k, n_dim = 1)
    G.fit(expr)
    [probs, clusters] = G.decode(expr)
    [probs, mixtures] = G.eval(expr)
    mean_as = argsort(G.means,0)
    for i in range(shape(mixtures)[0]):
        mixtures[i,:] = mixtures[i,squeeze(mean_as)]

    
    if draw:
        n = len(expr)
        xax = arange(n)[argsort(expr,0)]
        f = plt.figure(fig)
        f.clear()
        ax = f.add_axes([0,0,1,1])

        ct =mycolors.getct(k)
        cs, rs = [], []

        c2s  = []
        r2s = []
        x2s = []
        y2s = []
        for i in range(n):
            cs.append(ct[mean_as[clusters[i]]])
            rs.append(100)
            for j in range(k):
                mprob = mixtures[i,j]
                x2s.append(i)
                y2s.append(expr[i])
                c2s.append(ct[j])
                r2s.append(pow(exp(mprob),2)*100)
           
        x3s,y3s,c3s,r3s = [], [], [], []
        for i in range(n):
            probs = mixtures[i,:]
            cval = clusters[i]
            srt = argsort(probs)[::-1]
            maxval = probs[srt[0]]
            secval = probs[srt[1]]
            reliable = False
    
            if log(maxval) - log(secval) > cmin_diff: reliable = True

            x3s.append(i)
            y3s.append(expr[i])
            r3s.append(200)
            if not reliable:
                color = [0,0,0]
            else:
                color = ct[srt[0]]
            c3s.append(color)
    #ax.scatter(xax,expr,rs, color = cs)
        ax.scatter(x2s,y2s,r2s,edgecolor=c2s,facecolor = 'none')   
        ax.scatter(x3s,y3s,r3s,c3s)
    return mixtures

def expr_getonoff(expr_in):
    expr = array(expr_in)
    dev = std(expr)
    #expr = log(expr+dev)
    k = 3
    km = mlpy.Kmeans(k)
    n = len(expr)
    expr_2d = []
    for i in range(n):
        expr_2d.append(array([expr[i],0]))
    expr_2d = array(expr_2d)
    comp = km.compute(expr_2d)
    means = km.means
    compsort = arange(k)[argsort(map(lambda x: x[0],means))]

    n = len(expr)
    xax = argsort(expr)
    
    means = zeros(k)
    stds = zeros(k)
    for i in range(k):
        idxs = nonzero(equal(comp,i))[0]
        vals = array(expr)[idxs]
        means[i] = mean(vals)
        stds[i] = std(vals)
        

    f = plt.figure(1)
    f.clear()
    ax = f.add_axes([0,0,1,1])

    ct =mycolors.getct(k)
    cs, rs = [], []
    for i in range(n):
        cs.append(ct[compsort[comp[i]]])
        rs.append(100)
    ax.scatter(xax,expr,rs, color = cs)
    

    x0 = 0
    y0 = 0
    for i in range(k):
        ax.plot([x0,x0],[means[i] - stds[i], means[i]+ stds[i]]
                ,linewidth = 5, color = ct[compsort[i]])





#... default ordering for TFs
#right now, they are ordered simply by
#the order of their out degree
def tf_ordering(name = default_name, reset = 0 ):
    hardcopy = True
    try:
        if reset: raise Exception('compute')
        return nw.readnet(name,hardcopy = hardcopy)
    except Exception as e:
        if e.args[0] != 'compute': raise Exception()
        nw.claim_reset()

        #default keymap into tfspace
        keymap = na_tfkeys()
        aff = net_affinity()[0]
        deg_sort = np.argsort(np.sum(aff,0))[::-1]
        
        kordered = {}
        count = 0
        for d in deg_sort:
            k = keymap.keys()[keymap.values().index(d)]
            kordered[k] = count
            count+=1
        
        nw.writenet(name, kordered,hardcopy = hardcopy)
        return kordered



#... default ordering for genes
#right now, they are ordered simply by
#the order of their in degree
def gene_ordering(name = default_name, reset = 0 ):
    hardcopy =True
    try:
        if reset: raise Exception('compute')
        return nw.readnet(name,hardcopy = hardcopy)
    except Exception as e:
        if e.args[0] != 'compute': raise Exception()
        nw.claim_reset()

        #default keymap into tfspace
        keymap = na_gkeys()
        aff = net_affinity()[0]
        deg_sort = np.argsort(np.sum(aff,1))[::-1]
        
        kordered ={}
        count = 0
        for d in deg_sort:
            k = keymap.keys()[keymap.values().index(d)]
            kordered[k] = count
            count+=1

        #this way is faster but who cares?
        #klist = []
        #kvpairs = map(lambda x,y:[x,y], keymap.iteritems())
        #vals = map(lambda x: x[1],kvpairs)
        #vals_srt = argsort(vals)
        #keys = map(lambda x: x[0],kvpairs)
        #keys_srt = keys[vals_srt]
        
        
        
        nw.writenet(name, kordered,hardcopy = hardcopy)
        return kordered


