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

#Basically, getorset can do three things
#  1: SET the stored value for the output of a function
#  2: GET the stored value for the output of function
#  3: UPDATE a new value without calling the function



def _make_reset_level(**kwargs):
  '''Get the reset level for the current function'''
  reset = kwargs.get('reset', 0)
  if 'resets' in kwargs: 
    reset = kwargs['resets'].get( inspect.stack()[2][3], reset) 
    #raise Exception()
  return reset

def _make_sub_reset_level(**kwargs):
  reset = kwargs.get('reset', 0)
  if 'resets' in kwargs: 
    reset = kwargs['resets'].get( inspect.stack()[2][3], reset) 
  return mod(reset,2)

def sub_reconcile(kw_new, kw_old = {}, **kwargs):
  d0 = dict(kw_old)
  d0.update(kw_new)
  d0.update(kwargs)
  d0['reset'] = _make_sub_reset_level(**d0)
  return d0
sr = sub_reconcile 

def reconcile(kw_new, kw_old = {}, **kwargs):
  d0 = dict(kw_old)
  d0.update(kw_new)
  d0.update(kwargs)
  d0['reset'] = _make_reset_level(**d0)
  return d0
rc = reconcile
  
def getOrSet(function, **kwargs):
  reset = kwargs.get('reset', False)
  register = kwargs.get('register', 'a')
  name = kwargs.get('name',register)
  hardcopy = kwargs.get('hardcopy', True)
  np = kwargs.get('np', False)
  update = kwargs.get('update', None)
  on_fail = kwargs.get('on_fail', 'error')
  hard_reset = kwargs.get('hard_reset', False)

  caller_name = inspect.stack()[1][3]

  
  
  if update != None:
    out = update
    write(name = name, value = out,  hardcopy = hardcopy, np = np,
          register = register, caller_name = caller_name)
  elif not reset:
    out, sxs = read(name = name, hardcopy = hardcopy, np = np, 
                    register = register, caller_name = caller_name)
    if not sxs:
      if on_fail == 'compute': 
        print 'memo.py:\n  Fetch failed for {0}, name: {1}\n  "compute" flag is set'.\
            format(caller_name,name)
        reset = True
      else: assert 0, 'Data recovery failed for name ' + caller_name

  if reset:
    if hard_reset:
        user_inp = raw_input('''
This appears to be a hard function to compute ({0}:{1})
Really Reset? (y/n)
'''.format(caller_name, name))
        assert user_inp == 'y'

    out = function(**kwargs)
    write(name = name, value = out,  hardcopy = hardcopy, np = np,
          register = register, caller_name = caller_name)
  return out
    

def write(name = None,value = None, hardcopy = True, np = False, 
          register= 'a', caller_name = None):
    if name == None: name = register
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
def read(name = None, hardcopy = True,np = False, 
         register = 'a', caller_name = None):
    if name == None: name = register
    sxs = False
    out = -1

    if not caller_name:caller_name = inspect.stack()[1][3]
    try:
        lname = globals()['lastname_'+caller_name+register]
        if lname == name or name == register:
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
            try:
              f = open(path,'r')
            except Exception as e:
              return out, sxs
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
