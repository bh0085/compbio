import os
import numpy as np
import re
import itertools as it
import matplotlib.pyplot as plt
from scipy import sparse
from numpy import *
import pickle
import subprocess
import random
import network.netutils as nu

def SVM_predict(trg_n = 0, name = 'mRN', number = 0, option ='svd'):
    trgs, tfs = nu.parse_net(name,number)
    ts = nu.load_TS()
    cl = nu.load_CL()

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
    
