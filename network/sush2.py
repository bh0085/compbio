import os
import numpy as np
import re
import itertools as it
import matplotlib.pyplot as plt
from scipy import sparse
from numpy import *
import numpy as np
import pickle
import subprocess
import random
import network.netutils as nu
import colors
import mlpy
import mysvm

def get_membership(genes,method = 'svd',svd = 'U', n = 10):
    if method == 'identity':
        idenfun = lambda x,y: x==y and 1.0 or 0.0
        members = zeros((n,len(genes)))
        for i in range(n):
            members[i,i] = 1
        #members = array([[idenfun(i,j) for j in range(len(genes))] for i in range(n)])
    if method == 'svd':
        print 'getting membership values from the globals, "net_svd, net_svd_keys"'
        global net_svd
        global net_svd_keys
        
        if not net_svd: raise Exception('Sorry, looks like the global net_svd is unset')
        #columns of V give membership
        #columns of U give membership

        #grab svd eigenmatrices and transform them so that rows
        #give trg/tf membership vectors... (e.g that they are basis vectors)
        if svd == 'U':
            vecs = net_svd[0]
            vecs = vecs.T
            keys = net_svd_keys[0]
            
        else:
            vecs = net_svd[2]
            vecs = vecs.T
            keys = net_svd_keys[1]
        #each 'key' corresponds to a column in vecs
        #for each gene, find its position in the svd
        #fill the membership matrix with its membership in each tf.
        input_keys = map(lambda x: x[0],genes)    
        vec_column_indices = []
        ng = len(genes)
        for i in range(ng):
            vec_column_indices.append(  keys[input_keys[i]] )
        

        vecs = array(vecs[:,vec_column_indices])
        members = vecs[0:n,:]
        
    
    return members 
def get_trgs_tfs(name = 'fRN', number = 0):
    return nu.parse_net(name,number)

def get_chosen_xy(x_members, y_members, method = 'svd', choose = 'y', index = 0, n = 10, adj = None):
    '''Choose xy members for learning.
By default, assume we are working in the eigengene/eigentf space and 
specifying a gene of interest with choose = 'y' and index = 0.

N fixes the number of eigengenes that will be modeled.
Because the adjacency matrix is effectively the identity,
N also fixes the number of tfs that will be specified.
'''

    #choose a subset of x and y to learn on:
    #if method is identity, it is assumed that we can find a row (or rows)
    #in the 'choose' array with value C[index,j] = 1
    
    if choose == 'y':
        chooser = y_members
        nonchoo = x_members
        if adj != None: 
            a= adj
    else:
        chooser = x_members
        nonchoo = y_members
        if adj != None:
            a= adj.T

    #if choose == 'y' ('x'), a is a (transposed) version of the 
    #adjacency matrix as follows:

    #index i specifies a row of the matrix from which columns,
    #x (y) will be selected subject to the condition that column 
    #x (y) is nonzero in its 'i'th row.

    if method == 'svd':
        #in fact, this takes the same steps as the 'identity'
        #but the adjacency matrix is the identity. 
        max_repr = argsort( chooser[:,index] )[::-1][0:n]
        cset = max_repr
        nset = max_repr
    if method == 'identity':
        #find chooser members corresponding to index
        max_repr = nonzero( chooser[:,index] )[0]
        cset = max_repr

        #find elements adjacent to index
        adjacent = nonzero( a[index,:])[0]
        #identify all nonchooser members containing adjacent elements
        max_repr2 = nonzero( nonchoo[:,adjacent] )[0]
        nset = max_repr2
    else:
        print 'sorry no choice method besides svd and identity has been implemented yet.'
        
    if choose =='y':
        yset = cset
        xset = nset
    else:
        yset = nset
        xset = cset

    return xset,yset
    
    
    
    #returns a subset of x_members and a subset of y_members
def get_training_data(x_expr, y_expr, y = 0):
    #write_SVM takes training data in the format:
    #[...
    #[label_yi,[x1,x2,x3....]]
    #[...
    
    xset = x_expr[:,:]
    yset = y_expr[y,:]
    nx = shape(xset)[0]
    nt = shape(xset)[1]

    if len(shape(yset)) == 1:
        training_data =  [ [yset[i],xset[:,i]] for i in range(nt)]
    else:
        print 'multiclass learning not yet supported'
        raise Exception()
    return training_data, xset, yset



def run( trgs , tfs, name = 'fRN', number = 0, method ='svd',yindex = 0 ):
    
    trg_l = list(map(lambda x: [x[0],x[1]], trgs.items()))
    tf_l = list(map(lambda x: [x[0],x[1]], tfs.items()))
    ntf = len(tf_l)
    ntrg = len(trg_l)

    tf_keys = map(lambda x: x[0], tf_l)


    #choose a method: eigen or svd
    expression = 'time'
    trange = [0,100]

    #for embryonic expr only, use upper = 12
    time_upper = 30
    if expression == 'time':
        data =nu.load_TS()
        for k,v in data.items():
            #embryonic time points only
            data[k] = v[0:time_upper]
    else:
        data = nu.load_CL()

    #get a list of targets (things to be predicted)
    tf_trglist= []
    tf_idxlist= []

    trg_tflist = []
    trg_idxlist = []    

    trg_idxlist =range(len(trgs.keys()))
    for t in trg_idxlist:
        trg_tfs = []
        for tf in trg_l[t][1]['tfs']:
            tf_idx = tf_keys.index(tf)
            if not tf_idx in tf_idxlist:
                tf_idxlist.append(tf_idx)
                tf_trglist.append([])
                
            trg_tfs.append(tf_idx)
            tf_trglist[tf_idxlist.index(tf_idx)].append(t)
        trg_tflist.append(trg_tfs)

    adj = zeros((ntrg,ntf))
    for i in range(ntf):
        adj[tf_trglist[i][:],i] = 1


    #build membership lists for x (predictors)
    #and y (output) that will be used in learning.
    x_memberships = get_membership(tf_l, method = method, svd = 'V', n = len(tf_l))
    #get membership vectors according to the svd. search only for len(tf) vectors.
    y_memberships = get_membership(trg_l, method = method, svd = 'U', n = 100)

    tf_expr = array([ data[tfs.keys()[i]] for i in tf_idxlist])
    trg_expr = array([ data[trgs.keys()[i]] for i in trg_idxlist])
    nt = len(trg_expr[0])


    #get expression values x and y at the datapoints provided.
    y_expr =array([ [sum(trg_expr[:,i]*y_memberships[j,:]) for i in range(trg_expr.shape[1])] for j in range(y_memberships.shape[0]) ])
    x_expr =array([ [sum(tf_expr[:,i]*squeeze(x_memberships[j,:])) for i in range(tf_expr.shape[1])] for j in range(x_memberships.shape[0]) ])
    


    #look at a subset of xs and ys.
    xidxs,yidxs = get_chosen_xy(x_memberships, y_memberships, n = 10, index = yindex , method = method, adj = adj)
    x_expr = x_expr[xidxs,:]
    y_expr = y_expr[yidxs,:]
    
    draw_expr(x_expr)

    #for starters, predict each label seperately:
    ny = shape(y_expr)[0]
    i = 0
    while 1:
        if i >= ny: 
            i -=1
            break
        training_data, xt, yt = get_training_data(x_expr,y_expr, y = i)
        
        y = array(map(lambda x : x[0], training_data))
        x = array(map(lambda x:  x[1], training_data))
        actual = map(lambda x: x[0],training_data)


        #MLPY RIDGE REGRESSION
        regression = mlpy.RidgeRegression()
        regression.learn(x,y)
        reg_predictions = regression.pred(x)
        draw_prediction(reg_predictions,actual,fig=3)
        
        #MLPY ELASTIC PREDICTIONS
        #Note that in the paper, values of tau =2*10^-5, mu ~.005
        #work well.
        elastic = mlpy.ElasticNet(2e-5, 5e-3)
        elastic.learn(x,y)
        net_predictions = elastic.pred(x)
        draw_prediction(net_predictions,actual,fig=2)

        #SVM PREDICTIONS
        model = mysvm.SVM_train(training_data)
        predictions_file = mysvm.SVM_classify(training_data,model)
        predictions = mysvm.read_SVM(predictions_file)

        draw_prediction(predictions,actual,fig = 0)    
        draw_xy(x_expr,squeeze(y_expr[i,:]))
                                   

        r = raw_input('option?(p:prev,n:next,d:done)')
        if r =='p':
            i-=1
        elif r =='d':
            break
        else:
            i+=1
    

    return predictions, actual, x_expr, y_expr[y,:]


def draw_prediction(predictions, actual,fig = 0):
    f = plt.figure(fig)
    f.clear()
    ax = f.add_axes([0,0,1,1])
    xax = arange(0,len(predictions))
    ax.plot(xax,predictions)
    ax.plot(xax,actual)
    eps = std(actual)/2
    minline = actual - eps
    maxline = actual + eps
    ax.plot(xax,maxline,alpha = .3)
    ax.plot(xax,minline,alpha = .3)
    ax.fill_between(xax,predictions, maxline,
                    where = greater(predictions,maxline),
                    color = 'red',
                    interpolate = True)
    ax.fill_between(xax,predictions, minline,
                    where = less(predictions, minline),
                    color = 'blue',
                    interpolate = True)
    f.show()
def draw_expr(expr):
    f = plt.figure(5)
    f.clear()
    ax = f.add_axes([0,0,1,1])
    n = shape(expr)[0]
    ct = colors.getct(n)
    for i in range(n):
        e = expr[i]
        ax.plot(e,color = 'red')
def draw_xy(xset, yset):
    
    nx = shape(xset)[0]
    nt =shape(xset)[1]
    print nt, len(yset)

    ct = colors.getct(nx)

    f2 = plt.figure(1)
    f2.clear()
    ax2 = f2.add_axes([0,0,1,1])
    xs, ys, rs, cs = [], [], [], []
    for i in range(nx ):
        feature = xset[i]
        fmax = max(feature)
        for t in range(nt):
            xs.append(feature[t]/fmax)
            ys.append(yset[t])
            rs.append(20)
            cs.append(ct[i])

    ax2.scatter(xs,ys,rs,cs)
    
    f2.show()

    

