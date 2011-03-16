### orange_models.py --- 
from numpy import * 
import orange, orngSVM
import generic_model as gm

class OMultiMod(gm.GenericModel):
  def __init__(self):
    super(OMultiMod,self).__init__()
    self.regressors = []
    self.regressor_idxs = []
    self.domains = []
  def makeDomain(self,xsub):
    domain = orange.Domain([orange.FloatVariable() for x in range(len(xsub))],
                           orange.FloatVariable())
    return domain

  def trainingData(self,xsub, ysub,domain):
    data = orange.ExampleTable([orange.Example(domain,list(xs) + [0.0])
                                for xs in xsub.T]) 
    for d , y in zip(data,ysub): d.setclass(y)
    return data

  def testingData(self,xsub,domain):
    data = orange.ExampleTable([orange.Example(domain,list(xs) + [0.0])
                                for xs in xsub.T]) 
    return data

  def addRegressor(self,regressor, regressor_idxs, domain):
    self.regressors.append(regressor)
    self.regressor_idxs.append(regressor_idxs)
    self.domains.append(domain)

  def predict(self, xvals):
    ny = len(self.regressors)
    predictions = []
    for i in range(ny):
      xsubs = xvals[self.regressor_idxs[i]]   
      data = self.testingData(xsubs, self.domains[i])
      pred = array([self.regressors[i](e) for e in data])
      predictions.append(pred)
    return predictions
      


class KNNModel(OMultiMod):
  def __init__(self,params = []):
    super(KNNModel,self).__init__()
    self.k = ceil( float(params[0][0]) / params[0][1] * 10) if len(params) > 0 else 3

  def learn(self, xvals, yvals, couplings):
    ny = len(yvals)
    for i in range(ny):
      ysub = yvals[i]
      idxs = self.couplingFun(couplings, i)
      xsub = xvals[idxs]
      domain = self.makeDomain(xsub)
      data = self.trainingData(xsub, ysub, domain)
      learned_model = orange.kNNLearner(data, k = self.k)
      self.addRegressor(learned_model, idxs, domain)

class EpsSVMModel(OMultiMod):
  def __init__(self, params = [], pdict = {}):
    self.svmtype = orange.SVMLearner.Epsilon_SVR
    self.C = self.logParam(params, 0, 1., [-2,6])
    self.p = self.logParam(params, 1, 1.,  [-2,0])
    pdict['C'] =  self.C
    pdict['p'] = self.p
    pdict['pkeys'] = ['C', 'p']
    super(EpsSVMModel,self).__init__()
    
  def learn(self, xvals, yvals, couplings):
    ny = len(yvals)
    for i in range(ny):
      ysub = yvals[i]
      idxs = self.couplingFun(couplings, i)
      xsub = xvals[idxs]
      domain = self.makeDomain(xsub)
      data = self.trainingData(xsub, ysub,domain)
      l = orange.SVMLearner()
      l.svm_type = self.svmtype
      l.C = self.C
      l.p = self.p
      regressor = l(data)
      self.addRegressor(regressor, idxs, domain)

  

class NuSVMModel(OMultiMod):
  def __init__(self, params = [], pdict = {}):
    self.plist, self.pdict = \
        self.splitParams(params, 
                         pdict = pdict,
                         kw_idxs = dict(C = 0,
                                        nu = 1),
                         kw_bounds = dict(C=[-2.,2.],
                                          p = [-1., -.001]))
    self.svmtype =  orange.SVMLearner.Nu_SVR  
    super(NuSVMModel,self).__init__()
  
  def learn(self, xvals, yvals, couplings,params = {}, easy = True):
    ny = len(yvals)
    for i in range(ny):
      ysub = yvals[i]
      idxs = self.couplingFun(couplings, i)
      xsub = xvals[idxs]
      domain = self.makeDomain(xsub)
      data = self.trainingData(xsub, ysub,domain)
      if easy:
        l = orngSVM.SVMLearnerEasy(svmtype = self.svmtype)
      else:
        l = orngSVM.SVMLearner(svmtype = self.svmtype)
      #l.svm_type = self.svmtype
      l.C = self.pdict['C']
      l.nu = self.pdict['nu']
      regressor = l(data)
      self.addRegressor(regressor, idxs, domain)
  
