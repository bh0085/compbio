from numpy import *

#generic model doesn't do diddly yet.
class GenericModel(object):
  '''defines coupling and logspace helper functions that
  should be overriden by advanced classes'''

  def couplingFun(self,couplings, idx):
    '''
    Trivial coupling function couples and input and output variables
    whose coupling is nonzero.
'''
    return  nonzero(couplings[idx])[0]
  
  def logParam(self,params, idx, default,logrange = [-2,2]):
    '''
    Get a log transform of the given parameter at a specified index or
  return a default value.
  '''
    return logspace(logrange[0], logrange[1],params[idx][1])[params[idx][0]] \
        if len(params) >= idx + 1 else default
  
  def linParam(self,params, idx, default,linrange = [-2,2]):
    '''
    Get a log transform of the given parameter at a specified index or
  return a default value.
  '''
    return linspace(linrange[0], linrange[1],params[idx][1])[params[idx][0]] \
        if len(params) >= idx + 1 else default

  
  def splitParams(self,params, 
                  kw_idxs = {},
                  pos_bounds = [],
                  kw_bounds = {}, 
                  pos_logs = [], 
                  kw_logs = {},
                  pdict = {},
                  plist = []):
    '''
  Automatically split the parameters provided by the meta algorithm that is running parameter testing. 

accepts:
  params:     a list of tuples containing parameter values and maxima 
  pos_bounds: a list of positional pars with their max/minima
  kw_bounds:  a dict of kw pars with max/minima
  pos_logs:   space parameters logarithmically?
  kw_logs:    ''
  
returns:
  plist, pdict
  '''
    pdict['pkeys'] = kw_idxs.keys()
    idxs = range(len(params))
    for k,v in kw_idxs.iteritems():
      i = idxs.pop(idxs.index(v))
      bounds = kw_bounds.get(k, [-1.,1.])
      islog = kw_logs.get(k,True)
      if islog:
        pdict[k] = self.logParam(params,i, 1., bounds)
      else:
        pdict[k] = self.linParam(params,i, 1., bounds)
    
    for i in idxs:
      bounds = pos_bounds[i] if i < len(pos_bounds) else [-1.,1.]
      islog = pos_logs[i] if i < len(pos_logs) else True
      if islog:
        plist[i] = self.logParam(params, i, 1., bounds)
      else:
        plist[i] = self.linParam(params, i, 1., bounds)
    
    return plist, pdict
      
