#!/usr/bin/env python
'''
Utilities to read network data including patrick's unsup net 
and the timecourse/cell line data.

Reads stuff scattered around data/network into python variables.
Any function can take kwargs including 'reset'. Setting reset=1
will often fix problems with finding stuff that my utilies assume
to be precomputed.
'''

import compbio.utils.memo as mem
import compbio.config as config
from numpy import *
import numpy as np,  itertools as it, os, re


def getNet(**kwargs):
  '''Get the saved network from patrick's files.

  output: tuple of dicts keyed by gene/tf names

          trgs: {gname: {color:0.}{weights:[0....]}{tfs:['tfname']}
                 ...}
          tfs : {tfname:{color:0.}{weights:[0....]}{tgs:['tfname']}
                 ...}'''
  def setNet(**kwargs):
    fpath = config.dataPath('network/patrick/unsup_patrick.txt')
    TC = getTC( reset = mod(kwargs.get('reset',0),2))
    CL = getCL( reset = mod(kwargs.get('reset',0),2))
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
      if TC.has_key(tfs[i]) and CL.has_key(targets[i]):
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

    return  (trg_d, tf_d)
  return mem.getOrSet(setNet, **kwargs)
  pass

def getTC(**kwargs):
  '''Timecourse data

  output: dictionary keyed by gene names'''
  def setTC(**kwargs):
    f = open(config.dataPath('network/TC.geneexp')).read()
    elts =f.split('\n')
    seqdict = {}
    for e in elts:
        matches = list(re.finditer(re.compile('([^\s]+)'), e))
        if not len(matches): continue
        name = matches[0].group(1)
        seqdict[name] = []
        for i in matches[1:]:
            seqdict[name].append(float(i.group(1)))

    for k, v in seqdict.iteritems():
      seqdict[k] = array(v)
      
    return seqdict
  return mem.getOrSet(setTC, **kwargs)
  
def getCL(**kwargs):
  '''Cell line data
  
  output: dict keyed by gene names'''
  def setCL(**kwargs):
    f = open(config.dataPath('network/CL.geneexp')).read()
    elts =f.split('\n')
    seqdict = {}
    for e in elts:
        matches = list(re.finditer(re.compile('([^\s]+)'), e))
        if not len(matches): continue
        name = matches[0].group(1)
        seqdict[name] = []
        for i in matches[1:]:
            seqdict[name].append(float(i.group(1)))

    for k, v in seqdict.iteritems():
      seqdict[k] = array(v)
      
    return seqdict
  return mem.getOrSet(setCL, **kwargs)

def getSush(**kwargs):
  '''Get sushmita's regression weights and biases'''
  def setSush(**kwargs):
    path = config.dataPath('network/network_predmodel/regressionwts/fRN')
    bias_files = [ os.path.join( path, f) for f in os.listdir(path)  if 'bias' in f ]
    nw_files = [ os.path.join( path, f) for f in os.listdir(path)  if 'nw' in f ]
    
    bias_re = re.compile('(?P<gname>\S+)\s+(?P<level>\S+)')
    weight_re = re.compile('(?P<gname>\S+)\s+(?P<tfname>\S+)\s+(?P<level>\S+)')
    genes = {}
    for b in bias_files:
      for l in open(b).xreadlines():
	match = bias_re.search(l)
	genes[match.group('gname')] = dict(bias = match.group('level'))
    for n in nw_files:
      for l in open(n).xreadlines():
	match = weight_re.search(l)
	g = genes[match.group('gname')]
	g['tfs'] = g.get('tfs', []) + [match.group('tfname')]
	g['weights'] = g.get('weights', []) + [match.group('level')]
    return genes
  return mem.getOrSet(setSush, **kwargs)
    
		
