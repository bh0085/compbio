import os, sys, inspect
import utils.adapters as adapters
import mgr


model_names = [os.path.splitext(f)[0] 
               for f in os.listdir(os.path.join(\
      os.path.dirname(inspect.stack()[0][1]),
      'models/'))
          if os.path.splitext(f)[-1] == '.py' 
               and f[0:2] != '__' and f[0] != '.']

globals()['models'] =\
    __import__('models', globals(), locals(),\
                 fromlist=model_names)

di = [( k, models.__dict__[k] )
      for k in model_names]
globals().update(di)
globals()['getName'] = mgr.getName
