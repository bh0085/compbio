#!/usr/bin/env python 
from numpy import *
import numpy as np
import matplotlib.pyplot as plt
import compbio.learning.plots as lplots
import compbio.utils.plots as myplots
import compbio.utils.colors as mycolors
from regression_models import *
from orange_models import *
import itertools as it
import compbio.utils.pbar as pbar

class Learner():
  def __init__(self,x_data, y_data, coupling):
    '''Inputs: x_data, y_data: [Nx * Nt], [Ny * Nt] arrays'''

    assert shape(x_data)[1] == shape(y_data)[1]

    self.xvals = x_data
    self.yvals = y_data
    self.coupling = coupling
    self.nt = shape(x_data)[1]
    self.nx = shape(x_data)[0]
    self.ny = shape(y_data)[0]
    self.splitTraining()

  def splitTraining(self, seed= -1, train_frac = .6):
    if seed != -1: random.seed(seed)
    inds = arange(self.nt)
    random.shuffle(inds)
    self.train_idxs = inds[0:floor(train_frac *self.nt)]
    self.test_idxs = inds[floor(train_frac *self.nt):]
    
  def xyTrain(self):
    return self.xvals[:,self.train_idxs] ,\
        self.yvals[:,self.train_idxs]
  def xyTest(self):
    return self.xvals[:,self.test_idxs] ,\
        self.yvals[:,self.test_idxs]
  def setModel(self,model):
    self.model = model
  def learn(self):
    x,y = self.xyTrain()
    self.model.learn(x,y,self.coupling)
  def predictTraining(self):
    x, y = self.xyTrain()
    return self.model.predict(x)
  def predictTest(self):
    x,y = self.xyTest()
    return self.model.predict(x)
  def predictOther(self, x, y):
    return self.model.predict(x)
  def makePlots(self, name = 'No Name'):
    xtrain,ytrain = self.xyTrain()
    xtest,ytest = self.xyTest()
    ytrain_predicted = self.predictTraining()
    ytest_predicted = self.predictTest()
    

    ny = len(ytrain)
    f = plt.figure(1)
    f.clear()
    ax0 = f.add_subplot('211')  

    f1 = plt.figure(2)
    f1.clear()
    ax1 = f1.add_subplot('211')
    ct = mycolors.getct(ny)
    for actual,predicted, ax ,subtitle in [[ytest,ytest_predicted,ax0,'test predictions'],
                               [ytrain,ytrain_predicted,ax1,'training predictions']]:
      for i in range(len(actual)):
        lplots.plotPredictions(actual[i], 
                               predicted[i], 
                               ax,
                               color = ct[i]) 
        myplots.maketitle(ax, name, subtitle =subtitle )

  def testParams(self, model_class,
                 prediction= 'test', res = 10, dim = 1):

    #set up the grid of prediction parameters
    
    if len(shape(res)) == 0: res = (res,) * dim
    test_vals =list( it.product(*[[(x,r) for x in arange(r)] 
                                  for r in res]))

    #shall we predict holdout  or training set?
    if prediction == 'training':
      xyfun = self.xyTrain
      predictfun = self.predictTraining
    else:
      xyfun = self.xyTest
      predictfun = self.predictTest

    #set initial values for output variables
    ntest = len( xyfun()[1][0])
    rms = zeros(res)
    pdicts = reshape(array([{} for i in range(product(res))]),res)
    test_preds = reshape(\
      array([ zeros(ntest) for i in range(product(res))]),\
        concatenate([res+(ntest,)]))
    
    #test the learning method for each parameter
    bar = pbar.simple(len(test_vals))
    ct = 0
    for t in test_vals:
      ct += 1
      bar.update(ct)
      pdict = {}
      idxs = zip(map(lambda x: x[0], t))

      self.setModel(model_class(params = t,
                                pdict = pdicts[idxs][0]))
      self.learn()      
      xtest, ytest = xyfun()
      ypred = predictfun()
      rms[idxs] = std(ytest - ypred)
      test_preds[idxs] = ypred
    bar.finish()
      
    #create a dictionary of all of the output variables
    out = {}
    out['pdicts'] = pdicts
    out['test_rms'] = rms
    out['test_preds'] = test_preds
    out['actual_preds'] = ytest
    print ytest
    return  out
def randomDependent(nx, ny, nt):
  xvals = random.random((nx,nt))

  cxns = greater(random.random((ny,nx)),.2)
  cofs = cxns * random.random((ny,nx))
  yvals = np.dot(cofs, xvals)

  return xvals, yvals, cxns



def main():
  nx, ny, nt = 20,10,100
  xvals, yvals, couplings = randomDependent(nx,ny,nt)

  #xvals  = array([0,1,2,3,4])[newaxis,:]
  #yvals =array( [.1,1.3,2.5,2.9,4.1])[newaxis,:]
  #couplings  = [[1]]
  l = Learner(xvals,yvals,couplings)
  l.setModel(LassoRegressionModel)
  l.learn()
  l.makePlots(name = 'Test Data')

      
    
if __name__ == '__main__':
  main()
  exit(0)
