#Write variables to common blocks and the hardrive to
#increase efficiency of frequently called functions
import inspect
import os
import pickle
import cPickle
from numpy import *
import compbio.config as config

verbose_utils = True
def claim_reset():
    if verbose_utils:
        print '''Resetting:
   '''+inspect.stack()[1][3]+'''
'''

#if writenet is called without a name, use the value
#default.

def getOrSet(function, **kwargs):
  reset = kwargs.get('reset', False)
  name = kwargs.get('name','default')
  hardcopy = kwargs.get('hardcopy', True)
  register = kwargs.get('register', 'a')
  np = kwargs.get('np', False)
  caller_name = inspect.stack()[1][3]
  if not reset:
    out, sxs = read(name = name, hardcopy = hardcopy, np = np, 
                    register = register, caller_name = caller_name)
    assert sxs, 'Data recovery failed for name ' + caller_name
  else:
    out = function(**kwargs)
    write(name = name, value = out,  hardcopy = hardcopy, np = np,
          register = register, caller_name = caller_name)
  return out
    

def write(name = 'default',value = None, hardcopy = True, np = False, 
          register= 'a', caller_name = None):
    if not caller_name: caller_name = inspect.stack()[1][3]
    savename = caller_name + '_' +name + '.memo'
 
    globals()['lastname_'+caller_name +register] = name
    globals()['last_'+caller_name+register] = value

    if hardcopy:
        path = os.path.join(config.getTempPath(), savename)
        f = open(path,'w')
        if np:
            save( f,value)
        else:
            pickle.dump(value, f)
        f.close()

  

#if readnet is called without a name, it just
#returns the most recently globalized net.
def read(name = 'default', hardcopy = True,np = False, 
         register = 'a', caller_name = None):
    sxs = False
    out = -1

    if not caller_name:caller_name = inspect.stack()[1][3]
    try:
        lname = globals()['lastname_'+caller_name+register]
        if lname == name or name == 'default':
            return globals()['last_'+caller_name+register], True
        else: 
            raise Exception()
    except Exception as e:
        if not hardcopy:
            return out, sxs
        savename = caller_name +'_'+name + '.memo'
        path =  os.path.join(config.getTempPath(),savename)

        if hardcopy:
            path = os.path.join(config.getTempPath(), savename)
            f = open(path,'r')
            if np:
                out = load(f)
                sxs = True
            else:
                out = pickle.load( f)
                sxs = True
            f.close()

    
        globals()['lastname_'+caller_name +register] = name
        globals()['last_'+caller_name+register] = out

    return out, sxs
