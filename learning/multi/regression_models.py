import mlpy
from numpy import *

class LassoRegressionModel():
  def learn(self, xvals, yvals, couplings):
    ny = len(yvals)
    nx = len(xvals)
    self.regressions = []
    self.regression_idxs = []
    for j in range(ny):
      idxs = nonzero(couplings[j])[0]
      regression = mlpy.Lasso(25)
      regression.learn(xvals[idxs].T,yvals[j])
      self.regressions.append(regression)
      self.regression_idxs.append(idxs)
  def predict(self,xvals):
    ny = len(self.regressions)
    predictions = []
    for j in range(ny):
      idxs = self.regression_idxs[j]
      x_sub= xvals[idxs]
      predictions.append(self.regressions[j].pred(x_sub.T))
    return array(predictions)

class RidgeRegressionModel():
  def learn(self, xvals, yvals, couplings):
    ny = len(yvals)
    nx = len(xvals)
    self.regressions = []
    self.regression_idxs = []
    for j in range(ny):
      idxs = nonzero(couplings[j])[0]
      regression = mlpy.RidgeRegression(25)
      regression.learn(xvals[idxs].T,yvals[j])
      self.regressions.append(regression)
      self.regression_idxs.append(idxs)
  def predict(self,xvals):
    ny = len(self.regressions)
    predictions = []
    for j in range(ny):
      idxs = self.regression_idxs[j]
      x_sub= xvals[idxs]
      predictions.append(self.regressions[j].pred(x_sub.T))
    return array(predictions)
