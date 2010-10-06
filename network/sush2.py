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

def parse_net(name = 'mRN',number = 0):
    prefix = os.path.join('predmodel/regressionwts', name)
    nwdata = open(os.path.join(prefix,
                               'nw_'+str(number) + '.sif')).read()
    #A few functions defined here to be used later
    trgfun = lambda x: x[1]
    wtfun = lambda x:float( x[2] )
    tffun = lambda x: x[0]
    sigmafun = lambda x: 1 / (1 + np.exp(-x /1))

    r = re.compile('^[ ]*(?P<target>[^\s]+)\t(?P<tf>[^\s]+)\t(?P<weight>[^\s]+)'
                   ,re.M)
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


    #Extract a dictionary with information for each TF
    tf_d = {}
    for k, g in it.groupby(cat,key = lambda x: x[0]):
        l = list(g)
        tf_targets = map(lambda x: x[1],l)
        
        tf_d[k] = {'targets':tf_targets}

    return (trg_d, tf_d)


def SVM_predict(trg_n = 0, name = 'mRN', number = 0, option ='svd'):
    trgs, tfs = parse_net(name,number)
    ts = load_TS()
    cl = load_CL()

    nt = len(ts.items()[0][1])
    nc = len(cl.items()[0][1])

    #choose the first target and plot its expression
    k0 = trgs.keys()
    k = k0[trg_n]
    trg_tfs = trgs[k]['tfs']
    expr_median = mean(ts[k])

    
    keys , SVD = pickle.load(open('temp/svd.pickle'))
    U, S, V = SVD
    tkeys , tfkeys = keys
    ns = len(S)
    idx = tkeys[k]


    nvecs = 5
    randomize = False
    if not randomize:
        vecs = argsort(U[idx,:])
        vecs = vecs[::-1]
    else:
        vecs =range(len(S))
        random.shuffle(vecs)
    vecs = vecs[0:nvecs] 

    do_bias = True
    training_data = []
    tfseries = {}

    ts_data = []
    cl_data = []
    ts_actual = []
    cl_actual = []
    for i in range(nt):
        ts_data.append(map( lambda x: ts[x][i],tfkeys))
        ts_actual.append(ts[k][i])
    for i in range(nc):
        cl_data.append(map( lambda x: cl[x][i],tfkeys))
        cl_actual.append(cl[k][i])
    
    cl_actual = array(cl_actual)
    ts_actual = array(ts_actual)
    cl_data = array(cl_data)
    ts_data = array(ts_data)

    ntest = 5
    ntrain = nt - ntest
    shuffled = arange(nt)
    random.shuffle( shuffled)
    ts_tr_idxs = shuffled[0:ntrain]
    ts_te_idxs = shuffled[ntrain:]

    ts_training_actual = ts_actual[ts_tr_idxs]
    ts_testing_actual = ts_actual[ts_te_idxs]
    cl_testing_actual = cl_actual

    ts_training_data = array(ts_data)[ts_tr_idxs]
    ts_testing_data = array(ts_data)[ts_te_idxs]
    cl_data = array(cl_data)

    ts_labels = array(ts[k])
    cl_labels = array(cl[k])
    
    ts_training_labels = ts_labels[ts_tr_idxs]
    ts_testing_labels = ts_labels[ts_te_idxs]
    cl_testing_labels = cl_labels

    ts_train, tts = get_training_data(ts_training_data, ts_training_labels, array(V[:,vecs].T))
    ts_test, tes =get_training_data(ts_testing_data, ts_testing_labels, array(V[:,vecs].T))
    cl_test, ces =get_training_data(cl_data, cl_testing_labels, array(V[:,vecs].T))

    ##USE THIS CODE TO USE A LINEAR COMBINATION FROM THE SVD
    #
    #if option =='svd':
    #    for t in range(nt):
    #        label = ts[k][t]
    #        labeled = [label,[]]
    #        tf_timevals = map( lambda x: ts[x][t],tfkeys)
    #        for i in range(nvecs):
    #            tf_vec =array( squeeze(V[:,vecs[i]]))
    #    
    #            #For a sanity check, see where the current vectors net tfs lie in vec.
    #            ag =squeeze(( argsort(tf_vec)))
    #            #print ag
    #            for tf in trg_tfs:
    #                print str(float(nonzero(ag == tfkeys[tf])[0]))  + '%'
    #            dotproduct = tf_vec* tf_timevals
    #            labeled[1].append(sum(dotproduct))
    #            if not str(i) in tfseries.keys(): 
    #                tfseries[str(i)] = []
    #            tfseries[str(i)].append(sum(dotproduct))
    #        if do_bias:
    #            labeled[1].append(1)
    #        training_data.append(labeled)
    #else:
    #    #Extract the time series of each of the relevant tfs:
    #    
    #    for t in range(nt):
    #        label = ts[k][t]
    #        #if ts[k][t] > expr_median : label = 1 
    #        #else: label = -1
    #        labeled = [label, []]
    #        for tf in trg_tfs:
    #            labeled[1].append(ts[tf][t])
    #            if not tf in tfseries.keys(): tfseries[tf] = []
    #            tfseries[tf].append(ts[tf][t])
    #            if do_bias:
    #                labeled[1].append(1)
    #        training_data.append(labeled)


    model_file =SVM_train(ts_train)

    for i in range(3):
        if i == 0:
            dfile = ts_train
            profile = tts
            actual = ts_training_actual
        elif i ==1:
            dfile = ts_test
            profile = tes
            actual = ts_testing_actual
        elif i== 2:
            dfile = cl_test
            profile = ces
            actual = cl_testing_actual

        predictions_file = SVM_classify(dfile, model_file)
        predicted = open(predictions_file).read().split('\n')
        predicted = array(list(it.ifilter(lambda x: x!='',predicted)))
        predicted = array(map(lambda x: float(x), predicted))
        np = len(predicted)
    
        colors = []
        rs = []
        for j in range(np):
            c = [(predicted[j] - expr_median) *( actual[j] -expr_median) < 0,
                 0,
                 (predicted[j]  - expr_median)* (actual[j] - expr_median)>0]
            colors.append(c)
            rs.append(25)

    
        fig = plt.figure(i)
        plt.clf()
        
        xax = arange(np)
        yax = actual
        
        ax0 = plt.axes(frameon = False)
        ax0.plot(xax,yax,color = 'r')
        ax0.plot(xax,zeros(np)+expr_median,color = 'b')
        ax0.scatter(xax, predicted,rs,color = colors)

        for v in profile:
            ax0.plot(xax,v,color = 'w', alpha = .25, zorder = -1)

def get_training_data(data,labels,vectors, nv = 5,do_bias =True):
    #data should be keyed with vector names
    nt = len(data)
    tfseries = []
    training_data = []
    for t in range(nt):
        label = labels[t]
        labeled = [label,[]]
        this_data = squeeze(data[t,:])
        for i in range(nv):
            this_vec = squeeze(vectors[i,:])
            dotproduct = this_vec*this_data
            labeled[1].append(sum(dotproduct))

            #an extra variable to keep track of tf variation for plotting.
            if i >= len(tfseries): tfseries.append([])
            tfseries[i].append(sum(dotproduct))

        if do_bias:
            labeled[1].append(1)
        training_data.append(labeled)
    return training_data, tfseries

def SVM_train(data):
    training_file = write_SVM(data)
    model_file = 'temp/model.svm'
    subprocess.call('svm_learn -e .2 -c .2 -z r '+training_file+' '+model_file,shell = True)
    return model_file
    
def SVM_classify(data, model_file):
    unlabeled_data = list(data)
    for d in unlabeled_data: d[0] = 0
    datafile = write_SVM(unlabeled_data)
    predictions_file = 'temp/predictions.svm'
    subprocess.call('svm_classify ' + datafile + ' ' + model_file + ' ' + predictions_file
                    ,shell = True)
    return predictions_file

def write_SVM(data):
    #data should just be a list of size [t,(1,n)] containing
    #n features sampled in t instances with labels.
    fname = 'temp/svm.data'
    f = open(fname,'w')

    n = len(data[0][1])
    t = len(data)

    for e in data:
        line = str(e[0]) + ' '
        for i in range(n):
            line += str(i+1) + ':' + str(e[1][i]) + ' '
        f.write(line + '\n')
    f.close()

    return fname
    

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

def load_CL():
    f = pickle.load(open('CL.pickle'))
    return f
def load_TS():
    f = pickle.load(open('TC.pickle'))
    return f

def svd_net(name = 'mRN', number = 0):
    #possible names: bRN, kRN, fRN, mRN
    trgs, tfs = parse_net(name, number)
    print 'N Targets: ' + str(len(trgs.keys()))
    print 'N TFs: ' + str(len(tfs.keys()))
    
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
    m = len(kmap0.keys())
    #mmax = 500
    #m = mmax

    n = len(kmap1.keys())
    N = np.zeros([m,n])
    for k_trg,v_trg in trgs.items():
        i = kmap0[k_trg]
        for k_tf in v_trg['tfs']:            
            j = kmap1[k_tf]
            N[i][j] = 1.0

    import pickle
    svdname = name + '.svd'
    comp_svd = True
    if comp_svd:
        U, S, Vh = lin.svd(N)
        V = Vh.T
    else:
        U,S,V = pickle.load(open(svdname))

    #from mlabwrap import mlab
    #U2, S2, V2 = mlab.svd(N,nout = 3)
        
    V = mat(V)

    uvecs = U[:,:len(S)]
    vvecs = V
    dosave = ((kmap0,kmap1),(uvecs, S, vvecs))
    pickle.dump(dosave,open('temp/'+name+'_svd.pickle','w'))

    #prod = ( U2 * SARR) * (conj(V).T)
    ##prod = np.dot(U,np.dot(SARR,np.conj(V)))
    #
    #print len(nonzero(N)[0])
    #
    #f = plt.figure(0)
    #plt.clf()
    #ax = plt.axes(frameon = False)
    #
    ##f2 = plt.figure(1)
    ##ax1 = plt.axes(frameon = False)
    #
    #
    #print sum(S)
    #
    #plot_vals = True
    #if plot_vals:
    #    #ax.set_xlim([0,1])
    #    #ax.set_ylim([0,1])
    #    
    #    xax = np.arange(float(len(S)))/len(S)
    #    yax = S # / S.max()
    #    ax.plot(xax, yax)
    #else:
    #    #ax.imshow(U2)
    #    ax.imshow(prod)
    #   # print abs(N).mean()
    #   # print  abs(N-prod).mean()
    #    pass
    #    
    #print N
def svd_view():
    net_names = ['mRN','fRN','kRN','bRN']
    
    for i in range(len(net_names)):
        net_name = net_names[i]
        #if net_name != 'mRN': continue
        keys, svd = pickle.load(open('temp/'+net_name+'_svd.pickle'))
        U, S, V = svd
        kmap0 = keys[0]
        kmap1 = keys[1]
        m = len(kmap0.keys())
        n = len(kmap1.keys())
        print 'm: '+str(m) +  ' n: ' + str(n)
        
        nt = shape(U)[0]
        yvals = pow(reshape(mod(arange(m*n) , n),(m,n)),2)
        tots = sum( abs(U) *yvals,1)
        srt = argsort(tots)
        Usrt = U[srt,:]
        
        fig = plt.figure(i+4)
        fig.clear()
        ax = fig.add_axes([0,0,1,1],frameon = False,title = net_name)
        ax.imshow(abs(U[srt[0:m/2],:]),aspect = 'auto',interpolation = 'nearest')

        #ax.plot(tots)
