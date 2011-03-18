import compbio.config as config
import re
from numpy import *
import numpy as np
import compbio.utils.heatmap as hm
import numpy as np
import compbio.utils.memo as mem
import itertools as it
import matplotlib.pyplot as plt


def parseNet(num =  1,method = 'tree', reset = False):
  '''
  Get one of daniel's nets. Allowable numbers are 1-3 and allowable
  types are 'tree', 'svm'
'''

  def setNet(**kwargs):
    method =kwargs.get('method', 'tree')
    num = kwargs.get('num', 1)

    description_path = config.dataPath('::daniel/net%s_chip_features.tsv') % num
    data_path = config.dataPath('::daniel/informativeness/%s%s.txt') %(method,num)
    split_re = re.compile('\s')
    
    desc_open = open(description_path)
    description_cols = split_re.split(desc_open.readline().strip())
    description_vals = [split_re.split(l.strip()) for l in desc_open.readlines()]
    
    data_open = open(data_path)
    weight, tf, exp = zip(*[array(split_re.split(l.strip()), float) 
                           for l in data_open.readlines()])
    description = {}
    for i in range(len(description_cols)): 
      description[description_cols[i]] = [d[i] for d in description_vals]
      
    
    ntf = np.max(tf) + 1
    nexp = len(description.values()[0]) + 1

    raise Exception()
    
    grid = zeros((ntf,nexp))
    for vals in zip(weight,tf,exp): grid[vals[1], vals[2]] = float(vals[0])
    
    return grid
  return mem.getOrSet(setNet,
    reset = reset, 
    register = method,
    name = '%s%s' %(method,num),
    method = method,
    num = num)
  
  
def linearize(array):
  return reshape(array, product(shape(array)))

def compareall():
  f = plt.figure(0)
  plt.clf()
  for i in range(3):
    num = [1,3,4][i]
    
    ax = f.add_subplot('23%s'%(i+1))
    d1,corr1 = compare(num, ax = ax)
    ax2 = f.add_subplot('23%s'%(i+4))
    d2,corr2= compare(num,rnd_tfs = True, ax = ax2)
    ax.set_title('Net %s correlation svm vs. tree importance :\n  %s' % (num, corr1))
    ax2.set_title('Net %s (randomized) \n correlation: %s' % (num, corr2))
    ax.set_xlabel('tree')
    ax.set_ylabel('svm')
def compare(num = 1, reset = False, 
            ax = None, 
            rnd_tfs = False, rnd_exps  = False):

  if not ax:
    plt.clf()
    f = plt.figure(0)
    ax = f.add_subplot('111')
    
  
  net_svm = parseNet(num = num, method = 'svm', reset = reset)
  net_tree = parseNet(num = num, method = 'tree', reset = reset)

  if shape(net_svm) != shape(net_tree):

    print '''There is something weird about this data,
the shapes do not match'''
    print '...fixing?'
    if shape(net_svm) > shape(net_tree):
      net_tree2 = zeros(shape(net_svm))
      inds = list(it.product(*[arange(s) for s in shape(net_tree)]))
      net_tree2[zip(*inds)] =\
          net_tree[zip(*inds)]
      net_tree = net_tree2
    else:
      net_svm2 = zeros(shape(net_tree))
      inds = list(it.product(*[arange(s) for s in shape(net_svm)]))
      net_svm2[zip(*inds)] =\
          net_svm[zip(*inds)]
      net_svm = roll(net_svm2,0)

  if rnd_tfs:
    net_svm = net_svm[random.permutation(range(len(net_svm)))]
  if rnd_exps:
    net_svm = net_svm[:, random.permutation(range(shape(net_svm)[1]))]

    
  net_svm = linearize(net_svm)
  net_tree = linearize(net_tree)

  srted = argsort(reshape(net_tree,product(shape(net_tree))))
  srted = srted[nonzero(greater(net_svm[srted] * net_tree[srted], 0))[0]]
  
  #box_indices = numpy.dot(ndims**numpy.arange(ndims), binassign)


  print 'Correlation coefficient: '
  corr =    corrcoef(log(net_svm[srted]), log(net_tree[srted]))
  print corr
  ax.plot(log(net_svm[srted]), log(net_tree[srted]),
          marker = 'o',
          linestyle = '',
          color = 'red',
          alpha = .2)
  return [[net_svm, net_tree, srted],corr[0,1]]
  #hm.heatMap(grid, xlabel = 'TFs', ylabel = 'Chips')
  
  
  
