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
import netutils as nu
import compbio.utils.colors as colors
import mlpy
import compbio.learning.myrf as myrf
import compbio.learning.mysvm as mysvm
import compbio.utils.plots as myplots
import matplotlib.patches as patches
import compbio.network.s3 as s3
import netwriter as nw
import compbio.learning.gagd as gagd
import compbio.learning.synth_fann as sf

default_name = 'unsup'
min_tf_perc = .7

#defines a mapping from expression space
#to feature vector/prediction space...
#
#it may not wind up being useful -> for now
#use the membership method 'identity'.
#
#input: genes, a list of gene names to compute membership for
def get_membership(genes,method = 'identity',svd = 'U'):
    #membership vectors should have genes as entries... thus the output is [Nm, ng]
    if method == 'identity':
        members = zeros((len(genes), len(genes)))
        for i in range(len(members)): members[i,i] = 1
        #members = array([[idenfun(i,j) for j in range(len(genes))] for i in range(n)])
    #<<<NOTE: THE SVD MEMBERSHIP METHOD WILL NOT WORK AT THE MOMENT>>>
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

def get_chosen_xy(x_members, y_members, 
                  method = 'identity', choose = 'y', 
                  index = 0, adj = None):
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
        invec = zeros(shape(chooser)[1])
        invec[index]=1
        transformed_adj = dot(dot(nonchoo.T,adj.T),chooser)
        nonchoo_out = dot(transformed_adj, invec)
        
        ngreater = len(nonzero(greater(nonchoo_out,0))[0])    
        max_repr = argsort(nonchoo_out)[::-1][0:ngreater]
        cset = [index]
        nset = max_repr
    else:
        print 'sorry no choice method besides svd and identity has been implemented yet.'
        
    if choose =='y':
        yset = cset
        xset = nset
    else:
        yset = nset
        xset = cset

    return xset,yset
    

def get_subs(cluster = 0 ):
    trgss = get_trg_ss(cluster = cluster)
    tfss = get_tf_ss(  trgnames = trgss)
    return trgss, tfss


def tfl_hash():
    tgls,tfls = test_subs(nmax = 2000)
    
    tfl_hash = [list([]) for i in range(max(tfls)+1)]
    for i in range(len(tfls)):
        tfl_hash[tfls[i]].append(i)
    return tfl_hash
def test_subs(nmax = 2000,reset = 0, name = default_name):
    if not reset:
        out,sxs = nw.rn2(name)
        if not sxs: raise Exception()
        tgls,tfls = out
    else:
        tfls = []
        tgls = []
        for cidx in range(nmax):
            tgs, tfs = get_subs(cidx)
            tgls.append(len(tgs))
            tfls.append(len(tfs))
        nw.wn2(name, (tgls,tfls))
    return tgls, tfls

#get tfs correstponding to a cluster
def get_trg_ss(cluster = 0):
    #for now... fake it!
    trgs, tfs = nu.parse_net()
    #grab a list of 10 random trgs!

    sib_arr = s3.sib_lists()
    sibs = nonzero(sib_arr[cluster])[0]
    kidxs = nu.net_trg_keyidxs()
    trgl = list(array(kidxs.keys())[sibs])
    if len(trgl) == 0: raise Exception()
    return trgl

#get trgs corresponding to a cluster     
def get_tf_ss(cluster = 0, trgnames = None, basic = False):
    if basic:

        trgs, tfs = nu.parse_net()
    #grab a list of the tfs regulating 10 random trgs
        tfl = []
        for k in trgs.keys()[0:50]:
            item = trgs[k]
            tfs = item['tfs']
            for t in tfs:
                if not t in tfl:
                    tfl.append(t)
    else:
        min_regs = min_tf_perc
        nat = greater(nu.net_affinity(),0)
        tgkeys = nu.net_trg_keyidxs()
        tg_sub = nat[[tgkeys[k] for k in trgnames],:]
        mem_means = mean(array(tg_sub,float),0)
        tfkeys = nu.net_tf_keyidxs()
        tfl = []
        tfhash = ['']* len(tfkeys.keys())
        for k,v in tfkeys.items(): tfhash[v] = k
        for n in nonzero(greater(mem_means,min_regs))[0]:
            tfl.append(tfhash[n])
        
    return tfl

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





#get a bunch of clusters close to the specifications
def get_candidates(n_wanted = 40, ntf = 5, ntg = 10):
    tgls, tfls = test_subs()
    
    scores = power(array(tgls) - ntg,2) + power(array(tfls) - ntf,2)
    sarg = argsort(scores)

    tgkeys = nu.net_trg_keyidxs()


    #find matching clusters matching the requirements.
    found_hash = zeros(len(tgkeys.keys()))
    matches = []
    for i in sarg:
        if len(matches) > n_wanted: break
        tgs, tls = get_subs(i)
        m = np.mean(found_hash[[tgkeys[k] for k in tgs]])
        if m < .5:
            found_hash[[tgkeys[k] for k in tgs]] +=1
            matches.append(i)
        else:
            print m, 'redundant match: ', i
 
          
        
    return matches


def get_cluster_expr(cluster, random_tfs = False):
    trg_ssnames = get_trg_ss(cluster = cluster )
    tf_ssnames = get_tf_ss( trgnames = trg_ssnames)    

    if random_tfs:
        tf_kidxs = nu.net_tf_keyidxs()
        r =np.random.random_integers(0,len(tf_kidxs.keys()),len(tf_ssnames))
        tf_ssnames = []
        print 'Randomizing TFs'
        for i in r:
            tf_ssnames.append(tf_kidxs.keys()[i])

    n = nu.parse_net()
    ts = nu.load_TS()
        
    tf_vals = array([ts[k] for k in tf_ssnames]).T
    tg_vals = array([ts[k] for k in trg_ssnames]).T
    
    all_exprs = []
    for vstart in [tg_vals, tf_vals]:
        vals = vstart
        mvals = np.mean(vals,1)
        vals -= mvals[:,newaxis]
        svals = np.std(vals,1)
        if len(nonzero(equal(svals,0))[0]):
            raise Exception('only one tf...')
    
        vals /= svals[:,newaxis]
        
        for v in vals.T:
            v -= mean(v)
            v /= std(v)
        all_exprs.append(vals)
        

    #raise Exception()
    return (all_exprs, (trg_ssnames,tf_ssnames))

def non_normal_cluster_expr(trg_ssnames, tf_ssnames, ctype = False, random_tfs = False):
    if random_tfs:
        tf_kidxs = nu.net_tf_keyidxs()
        r =np.random.random_integers(0,len(tf_kidxs.keys()),len(tf_ssnames))
        tf_ssnames = []
        print 'Randomizing TFs'
        for i in r:
            tf_ssnames.append(tf_kidxs.keys()[i])

    n = nu.parse_net()
    ts = nu.load_TS()
    if ctype:
        cl = nu.load_CL()
        tf_vals = array([ts[k] + cl[k] for k in tf_ssnames]).T
        tg_vals = array([ts[k] + cl[k] for k in trg_ssnames]).T
    else:
        tf_vals = array([ts[k]  for k in tf_ssnames]).T
        tg_vals = array([ts[k]  for k in trg_ssnames]).T        
    
    tf_vals -= np.mean(tf_vals,0)[:]
    tg_vals -= np.mean(tg_vals,0)[:]
    tf_vals /= np.std(tf_vals,0)[:]
    tg_vals /= np.std(tg_vals,0)[:]
        
    return [tg_vals, tf_vals]

def normalize_cluster_expr(trg_ssnames, tf_ssnames, ctype = False, random_tfs = False):
    if random_tfs:
        tf_kidxs = nu.net_tf_keyidxs()
        r =np.random.random_integers(0,len(tf_kidxs.keys()),len(tf_ssnames))
        tf_ssnames = []
        print 'Randomizing TFs'
        for i in r:
            tf_ssnames.append(tf_kidxs.keys()[i])

    n = nu.parse_net()
    ts = nu.load_TS()
    if ctype:
        cl = nu.load_CL()
        tf_vals = array([ts[k] + cl[k] for k in tf_ssnames]).T
        tg_vals = array([ts[k] + cl[k] for k in trg_ssnames]).T
    else:
        tf_vals = array([ts[k]  for k in tf_ssnames]).T
        tg_vals = array([ts[k]  for k in trg_ssnames]).T        

    all_exprs = []
    for vstart in [tg_vals, tf_vals]:
        vals = vstart
        mvals = np.mean(vals,1)
        vals -= mvals[:,newaxis]
        svals = np.std(vals,1)
        if len(nonzero(equal(svals,0))[0]):
            raise Exception('only one tf...')
    
        vals /= svals[:,newaxis]
        
        for v in vals.T:
            v -= mean(v)
            v /= std(v)
        all_exprs.append(vals)
        

    #raise Exception()
    return all_exprs
    

def viewclusters(cands, fig = 5):
    
    #clusters = get_candidates(get_candidates(10))
    
    f = plt.figure(fig)
    f.clear()
    ax1 = f.add_subplot(111)


    f2 = plt.figure(fig+1)
    f2.clear()
    ax2 = f2.add_subplot(111)

    for index in cands:
        trg_ssnames = get_trg_ss(cluster = index )
        tf_ssnames = get_tf_ss( trgnames = trg_ssnames)    

        n = nu.parse_net()
        ts = nu.load_TS()
        
        tf_vals = array([ts[k] for k in tf_ssnames]).T
        tg_vals = array([ts[k] for k in trg_ssnames]).T
        

        vals = tf_vals
        mvals = np.mean(vals,1)
        vals -= mvals[:,newaxis]
        svals = np.std(vals,1)
        if len(nonzero(equal(svals,0))[0]):
            raise Exception('only one tf...')

        vals /= svals[:,newaxis]

        for v in vals.T:
            v -= mean(v)
            v /= std(v)
        

        ax1.plot(vals)
        
        break
        
        
        #tf_vals = array([ts[k] for k in tf_ssnames]).T
        #x2.plot(np.std(tf_vals,1))
        
      
def run(  method ='identity',index = 0, reset = 0, 
          nxmax = 100 , 
          binary_x = False, binary_y = False, 
          expression = 'time' ,
          cluster_idx = 0,
          lrn = 'tree',
          showall = False,
          tgonly = False,
          randomize_tfs = False,
          ctfs = 5,
          ctgs = 5,
          cofs = 1,
          do_normalize_cluster = True,
          cluster_tfs = True,
          verbose_expr_labels = False,
          ctype = False):
    '''
sush2.run:

run a selected learning algorithm for  a cluster.

KEYWORDS:

index  [0]: select a tf/target to model from the cluster
method ['identity']: a membership method
multi  [False]: meaningless
nxmax  [3]: max cluster members
binary_x: model x data as binary
binary_y: model y data as binary
expression ['time']: which expression series to use
cluster_idx: not yet implemented

reset

'''

    #Data assembly:
    #
    #1: Grab a list of genes of interest and 
    #   corresponding expression vectors
    #
    trg_kidxs = nu.net_trg_keyidxs()
    tf_kidxs = nu.net_tf_keyidxs()
    #
    #retrieve the list of trg/tf names present in a given cluster.
    #note that at the moment, these are fake functions that just give back
    #a little list of trgs and all of their associated TFs
    #
    #--CLUSTERS USED--


    cands = get_candidates(10,ctfs,ctgs)
    cidx = cands[cofs]
    trg_ssnames = get_trg_ss(cluster = cidx )
    tf_ssnames = get_tf_ss(cluster = cidx , trgnames = trg_ssnames)
            
    if cluster_tfs:
        tf_ssnames = get_tf_ss(cluster = cidx , trgnames = trg_ssnames)
    else:
        tgs, tfs = nu.parse_net()
        tg_specific = trg_ssnames[cluster_idx]
        trg_tfs = tgs[tg_specific]['tfs']
        tf_ssnames = trg_tfs


    if randomize_tfs:
        r =np.random.random_integers(0,len(tf_kidxs.keys()),len(tf_ssnames))
        tf_ssnames = []
        print 'Randomizing TFs'
        for i in r:
            tf_ssnames.append(tf_kidxs.keys()[i])

    trg_ssidxs = array([trg_kidxs[name] for name in trg_ssnames])
    tf_ssidxs = array([tf_kidxs[name] for name in tf_ssnames])
    #
    #2: Project expression data onto membership vectors
    #
    #--EXPR CLUSTERING--
    #4: Grab a list of 'membership vectors' which
    #   translate genes to x and y in the machine learning problem
    #   data merging has not yet been implemented but should be quite simple
    #
    x_memberships = get_membership(tf_ssnames, method = method)
    y_memberships = get_membership(trg_ssnames, method = method)



    if do_normalize_cluster:
        exprtype = 'clustered'
    else:
        exprtype = 'standard'

    if exprtype == 'standard':
        all_expr = non_normal_cluster_expr(trg_ssnames, tf_ssnames,ctype = ctype)
    else:
        all_expr = normalize_cluster_expr(trg_ssnames, tf_ssnames,ctype = ctype)
        
    tg_expr, tf_expr = all_expr
    x_expr = array((tf_expr)).T
    y_expr = array((tg_expr)).T


    show_clustered_expr(y_expr,x_expr, trg_ssnames, tf_ssnames,fig = 8)    

    nx, npertg = shape(x_expr)
    x_all, y_all = fold_expr(x_expr, y_expr)
    nx, nt_folded = shape(x_all)
    train_idxs, test_idxs = [],[]

    nt = npertg
    if ctype:
        nt -= 4
    tginds = range(cluster_idx *npertg,(cluster_idx*npertg)+npertg)
    
    cinds = []
    for i in range(nt_folded):

        if (divmod(i,npertg))[1] >= npertg - 4:
            cinds.append(i)

    for i in range(nt_folded):
        if ctype:
            if i in cinds and i in tginds:
                test_idxs.append(i)
        else:
            if i in tginds[:-4]:
                test_idxs.append(i)
        if tgonly:
            if i in tginds[:-4]:
                train_idxs.append(i)
        else:
            if not (i in tginds) and not (i in cinds):
                train_idxs.append(i)
        


    print 'N_TRAIN' , len(train_idxs)
    expr_fig = 0
    draw_expr(x_expr, y_expr, expr_fig = expr_fig)

    if lrn =='svm':
        model = learn_svm( x_all, y_all,
                           train_idxs = train_idxs,
                           test_idxs = test_idxs,
                           binary_x = binary_x,
                           binary_y = binary_y)
        predictions = run_svm((x_all.T)[test_idxs].T , y_all[test_idxs], model)
    if lrn in ['knn','tree','forest']:

        #pred = myrf.run_tree(x_all,y_all, train_idxs, test_idxs)
        #raise Exception()

        all_ex = myrf.get_ex(x_all,y_all)
        train_ex = all_ex.getitems([int(x) for x in train_idxs])    
        test_ex  = all_ex.getitems([int(x) for x in test_idxs])    

        #test_ex = myrf.examples_from_inds(x_all,y_all,test_idxs)
        #cl_ex = myrf.examples_from_inds(x_all,y_all,cl_idxs)
        model = myrf.OLearn(lrn, train_ex, test_ex = test_ex)
        predictions = model.predictions(test_ex)

    if lrn == 'nn':

        nhc = 2
        ntg = 2
        ntf_s = 2
        max_tfu = 2
        gf = sf.genfann(nhc,ntg,ntf_s, [ max_tfu for i in range(ntg) ] )
        xs, ys = sf.synth_data(ntg,max_tfu,ntf_s)
        g, ga = gf.sample_genome()
        gf.init_net()
        gf.make_cxns_from_genome(g)
        #gf.net_from_cxns(hidden_cxns,output_cxns)

        net = gf.mynn.net
        
        f = plt.figure(0)
        f.clear()
        ax = f.add_subplot(121)
        myplots.draw_pb(ax,net)
        myplots.hideaxes(ax)
        myplots.maketitle(ax,'GANN')
        
        gf.set_data(xs.T,ys.T)
        gf.set_trainer()
        gf.train()


        ax2 = f.add_subplot(122)
        myplots.draw_pb(ax2,net)
        myplots.hideaxes(ax2)
        myplots.maketitle(ax2,'GANN')


        


        return
        raise Exception()

 
        


        
        raise Exception()

        #igrps = [ arange(2)+2*i for i in range(3) ]
        #igrps = [ 
        
        raise Exception()
        gf.train()

        raise Exception()
        #gagd.MyFANN(x_all.T,y_all[newaxis,:].T,train_idxs)

    actual = y_all[test_idxs]
    
    showall = True
    if showall:
        if verbose_expr_labels:
            names = tf_ssnames
        else:
            names = None
        draw_svm(x_all[:,test_idxs],actual, predictions, f = expr_fig,names = names)

    print predictions
    print actual

    if ctype:
        forstring = 'CL Data'
    else:
        forstring = 'TS Data'
        
    namestr = trg_ssnames[cluster_idx]
    subt = 'TFs: '+','.join(tf_ssnames)

    if randomize_tfs:
        title = 'Random TF Predictions ' + forstring + ', ' +namestr
        fnum = 5
    else:
        if cluster_tfs:
            title = 'Network Cluster TF Predictions'+ forstring + ', ' +namestr
        else:
            title = 'Network UnClustered TF Predictions'+ forstring + ', ' +namestr
            
        fnum = 6

    msecov = draw_prediction(predictions,actual,fig=fnum, 
                    title = title,
                    subt = ','.join(tf_ssnames))  

    print msecov
    return msecov

def show_clustered_expr(tge,tfe,tgnames, tfnames, nrml = True,fig = 8):
    f1 = plt.figure(fig)
    f2 = plt.figure(fig + 1)
    f1.clear()
    f2.clear()

    ax1 = f1.add_subplot(111)
    ax2 = f2.add_subplot(111)
    
    
    tgct = colors.getct(len(tgnames))
    tfct = colors.getct(len(tfnames))
    for i in range(len(tge)):
        ax1.plot(tge[i],color = tgct[i])
    myplots.color_legend(f1,tgct,tgnames, ax = ax1,pos = 4)
    tstr = 'Target Expression Levels' 
    if nrml: tstr += '(Normalized)'
    myplots.maketitle(ax1,tstr)

    for i in range(len(tfe)):
        ax2.plot(tfe[i],color = tfct[i])
    myplots.color_legend(f2,tfct,tfnames, ax = ax2,pos = 4)
    tstr = 'TF Expression Levels' 
    if nrml: tstr += '(Normalized)'
    myplots.maketitle(ax2,tstr)

    #raise Exception()

def normalize(x, binary = False):
    xdat = array(x)
    if binary:
        for elt in x:
            xdat[less_equal(x,0)] = 0
            xdat[greater(x,0)] = 1
    else:
        for elt in xdat:
            elt -= median(elt)
            elt /= std(elt)
    return xdat


def fold_expr(x_expr, y_expr):
    #if method == foldy, treat all ys as the same label and concatenate x.
    method = 'foldy'
    yout = []
    xout = []
    if method =='foldy':
        for y in y_expr:
            xout.extend(x_expr.T)
            yout.extend(y)


    xout  = array(xout)
    yout = array(yout)
    return xout.T, yout
    


def get_projected_exprs(members, exprs,binary = False,  method = 'linear'):
    '''currently only linear and median methods implemented
median is automatically toggled for binary data

'''
    if binary:
        method = 'median'
    
    if method =='linear':
        xpr = np.sum(members[:,:,newaxis] * exprs[newaxis,:,:],1)
    elif method =='median':
        xpr = np.median(members[:,:,newaxis] * exprs[newaxis,:,:],1)
    else:
        raise Exception('method not yet implemented')

    return xpr
    


def run_regression(x,y):
        #MLPY RIDGE REGRESSION
        regression = mlpy.RidgeRegression()
        regression.learn(x,y)
        reg_predictions = regression.pred(x)
        draw_prediction(reg_predictions,y,fig=3,subt = subt)

def run_elastic(x,y): 
        #MLPY ELASTIC PREDICTIONS
        #Note that in the paper, values of tau =2*10^-5, mu ~.005
        #work well.
        elastic = mlpy.ElasticNet(2e-5, 5e-3)
        elastic.learn(x,y)
        net_predictions = elastic.pred(x)
        draw_prediction(net_predictions,actual,fig=2)

def learn_svm(x, y, train_idxs = None, test_idxs = None,
              binary_y = False, 
              binary_x = False, 
              expr_figa =0, 
              expr_figb =0):

    tdat = ((x.T)[train_idxs,:],y[train_idxs])     
    pstring = mysvm.SVM_pstring(e = .2, binary = True, degree = 2, kernel = 'polynomial')
    model = mysvm.SVM_train(tdat, pstring = pstring)
    return model

def run_svm(x_expr,y_expr, model):
    tdat = (x_expr.T,squeeze( y_expr.T))
    predictions_file = mysvm.SVM_classify(tdat,model)
    predictions = mysvm.read_SVM(predictions_file)
    return predictions

def draw_svm(x_expr, predictions,actual, f = 0, names = None):

    fig = plt.figure(f)
    nx = len(x_expr)

    thr = .05
    pcols = map(lambda x : x <= -thr and 'blue' or \
                    x >= thr and 'red' or \
                    'black',
                predictions)    
    ycols = map(lambda x : x <= -thr and 'blue' or \
                    x >= thr and 'red' or \
                    'black',
                actual)
    mcols = map(lambda x, y: y == x and 'none' or \
                    'black', pcols, ycols)

    pcount = 0
    pmax = (power(nx,2) - nx) /2
    for i in range(nx):
        for j in range(nx):
            ax_r = [float(i)/nx, float(j)/nx,float(1)/nx,float(1)/nx]
            ax = fig.add_axes(ax_r,frameon = False)
            
            rsml = 20
            rbig = rsml*2.5

            ax.scatter( x_expr[i],
                x_expr[j],
                rsml,
                color = pcols
                )
            ax.scatter(x_expr[i],
                       x_expr[j],
                       rbig, 
                       edgecolor = ycols,
                       color = 'none')
            scatter_errors = False
            if scatter_errors:
                ax.scatter(x_expr[i],
                           x_expr[j],
                           rbig*2, 
                           color = mcols,
                           zorder = -100,
                           edgecolor = 'none')
            myplots.hideaxes(ax)
            if names:
                namearr = [names[i],names[j]]
            else:
                namearr= [str(i),str(j)]
            if names: alpha = 1.0
            else: alpha = .6
            myplots.maketitle(ax,' vs '.join(namearr), alpha = alpha)
                
def draw_prediction(predictions, actual,fig = 0, match_mean = True,title = '',subt = ''):

    
    f = plt.figure(fig)
    f.clear()
    ax = f.add_subplot(111)
    xax = arange(0,len(predictions))
    
    p2 = predictions - mean(predictions) 
    if std(p2) != 0: p2 = p2/std(p2) *std(actual)
    p2 = p2 + float(np.mean(actual))
    ax.plot(xax,p2)

    mse_zscore =np.sum( power(( p2 - actual),2))
    cov =np.sum( np.corrcoef( p2, actual)[0,1])
    if cov != cov: cov = 0

    ax.plot(xax,actual)
    eps = std(actual)
    minline = actual - eps
    maxline = actual + eps
    ax.plot(xax,maxline,alpha = .3)
    ax.plot(xax,minline,alpha = .3)
    
    ax.fill_between(xax,p2, maxline,
                    where = greater(p2,maxline),
                    color = 'red',
                    interpolate = True)
    ax.fill_between(xax,p2, minline,
                    where = less(p2, minline),
                    color = 'blue',
                    interpolate = True)
    myplots.maketitle(ax,title, subtitle = 'Validation MSE: '+str(round(mse_zscore,3))+'\nValidation Correlation: '+str(round(cov,3)))
    myplots.label_lr(ax,subt)
    f.show()
    return [mse_zscore,cov]
def draw_expr(x_expr,y_expr, expr_fig = 0 ):

    f = plt.figure(expr_fig)
    f.clear()
    nx = len(x_expr)
    cc = np.corrcoef(x_expr)
    #cc = (cc - np.min(cc)) / (np.max(cc) - np.min(cc))
    #colors = array([1,1,1])[newaxis,newaxis,:] + array([0,-1,-1])[newaxis,newaxis,:]*cc[:,:,newaxis]
    colors = color_array(cc)
    myplots.scatter_backdrop(f,colors)
    
def draw_xy(xset, yset):
    
    nx = shape(xset)[0]
    nt =shape(xset)[1]
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

    


def color_array(arr, maxabs = 1):
    maxabs = np.max(np.abs(arr))
    a2 = array(arr)
    a2/=maxabs

    slf = .75
    a2great = array([-1.0,-1.0,0.0])* ((a2* greater(a2,0      ))[:,:,newaxis])*slf
    a2less  = array([0.0,-1.0,-1.0])* ((-1*a2* less_equal(a2,0))[:,:,newaxis])*slf
    colors = array([1,1,1])[newaxis,newaxis,:] + a2great + a2less
    
    return colors
