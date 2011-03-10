### orange_models.py --- 
from numpy import * 
import orange

class KNNModel(k = 4):
  def learn(self, xvals, yvals, couplings):
    self.regressors = []
    self.regressor_idxs = []
    self.domains = []
    ny = len(yvals)

    for i in range(ny):
      ysub = yvals[i]
      idxs = nonzero(couplings[i])[0]
      xsub = xvals[idxs]
      domain = orange.Domain([orange.FloatVariable() for x in range(len(xvals))],
                           orange.FloatVariable())
      data = orange.ExampleTable([orange.Example(domain,list(xs) + [0.0])
                                  for xs in xsub.T]) 
      for d , y in zip(data,ysub): d.setclass(y)
      self.regressor_idxs.append(idxs)
      self.regressors.append(orange.kNNLearner(data, k = k))
      self.domains.append(domain)


  def predict(self,xvals):
    ny = len(self.regressors)
    predictions = []

    for i in range(ny):
      xsubs = xvals[self.regressor_idxs[i]]      
      data =   orange.ExampleTable([orange.Example(self.domains[i],list(xs) + [0.0])
                                    for xs in xsubs.T]) 
      pred = array([self.regressors[i](e) for e in data])
      predictions.append(pred)

    return predictions
      
