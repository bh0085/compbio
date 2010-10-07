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
        
        tf_d[k] = {'targets':map(trgfun,l),
                   'weights':map(wtfun,l)}

    return (trg_d, tf_d)


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
