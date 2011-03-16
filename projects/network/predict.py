#!/usr/bin/env python

import compbio.learning.multi.learner as l
import compbio.projects.network.io as netio
from numpy import *
import matplotlib.pyplot as plt
import compbio.learning.multi.orange_models as om
import compbio.learning.multi.regression_models as rm
import compbio.utils.colors as mycolors
import compbio.utils.plots as myplots
import compbio.utils.heatmap as hm

def learnGene(gene_name = 'FBgn0014931', 
              model = None ):
  if model == None: model = l.KNNModel()
  xvals,yvals,couplings = gVals(gene_name)
  learner = l.Learner(xvals, yvals, coupling)
  learner.setModel(model)
  learner.learn()
  learner.makePlots(name = gene_name)

def gVals(gene_name):
  TC = netio.getTC()
  trgs, tfs = netio.getNet()
  g = trgs[gene_name]
  gtfs = g['tfs']
  yvals = array([TC[gene_name]])
  xvals = array([TC[f] for f in gtfs])

  print 'TG Matrix shape: {0}'.format(shape(xvals))
  print 'TF Matrix shape: {0}'.format(shape(yvals))

  coupling = ones(len(xvals))[newaxis,:]
  return xvals, yvals, coupling
  

def heatMapGene(gene_name = 'FBgn0014931', 
                model_class = None,
                res = 5,
                prediction ='training'):
  plt.clf()
  if model_class == None: model_class = om.NuSVMModel
  xvals,yvals,coupling = gVals(gene_name)
  learner = l.Learner(xvals,yvals,coupling)
  vals = learner.testParams(model_class, 
                            prediction=prediction
                            ,res = res, dim = 2)
  err = vals['test_rms']
  annotations = vals['pdicts']
  f=plt.gcf()
  ax = f.add_subplot('211')

  ax2 = f.add_subplot('212')
  ax = hm.heatMap(err, annotations,axes = ax)
  
  myplots.maketitle(ax,  'gene: {0}'.format(gene_name), 
                    'heatmap for different learning parameters')
  
  preds = vals['test_preds']
  best_p = preds[unravel_index(argmin(vals['test_rms']),
                               shape(preds)[:2])]
  worst_p = preds[unravel_index(argmax(vals['test_rms']),
                               shape(preds)[:2])]
  ax2.plot(worst_p,
           linestyle = ':',
           linewidth = 4 ,
           color = 'blue')
  ax2.plot(best_p,
           linestyle = ':',
           linewidth = 4,
           color = 'red')
  ax2.plot(vals['actual_preds'][0])

#ax2.plot(sorted(reshape(err, (prod(shape(err))))))
  
  #plt.imshow(err, interpolation = 'nearest', cmap = mycolors.blackbody())
  
def cycleGenes(rng = [0,10], model_class = None):
  trgs, tfs = netio.getNet()
  for k in trgs.keys()[rng[0]:rng[1]]:
    print 'gene: {0}'.format(k)
    if model_class:
      model = model_class()
    else: model = None
    learnGene(k, model = model)

    f1 = plt.figure(1)
    plt.draw()
    f1.canvas.draw()
    f2 = plt.figure(2)   
    plt.draw()
    f2.canvas.draw()

    a = raw_input('    next gene?: (y/n) [y]    ')

    if a == 'n':
      break
