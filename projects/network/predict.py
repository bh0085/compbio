#!/usr/bin/env python

import compbio.learning.multi.learner as l
import compbio.projects.network.io as netio
from numpy import *
import matplotlib.pyplot as plt

def learnGene(gene_name = 'FBgn0014931', 
              model = None ):
  if model == None: model = l.KNNModel()

  TC = netio.getTC()
  trgs, tfs = netio.getNet()
  g = trgs[gene_name]
  gtfs = g['tfs']
  yvals = array([TC[gene_name]])
  xvals = array([TC[f] for f in gtfs])

  print 'TG Matrix shape: {0}'.format(shape(xvals))
  print 'TF Matrix shape: {0}'.format(shape(yvals))

  coupling = ones(len(xvals))[newaxis,:]
  learner = l.Learner(xvals, yvals, coupling)
  learner.setModel(model)
  learner.learn()
  learner.makePlots(name = gene_name)
  
def cycleGenes(rng = [0,10]):
  trgs, tfs = netio.getNet()
  for k in trgs.keys()[rng[0]:rng[1]]:
    print 'gene: {0}'.format(k)
    learnGene(k)

    f1 = plt.figure(1)
    plt.draw()
    f1.canvas.draw()
    f2 = plt.figure(2)   
    plt.draw()
    f2.canvas.draw()

    a = raw_input('    next gene?: (y/n) [y]    ')

    if a == 'n':
      break


def main():
  model = l.LinearRegressionModel()
  predictWithMethod(model)

if __name__ == '__main__':
  main()
  exit(0)
