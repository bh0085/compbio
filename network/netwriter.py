#Write nets to common blocks and the hardrive to
#increase efficiency of frequently called functions
import inspect
import os
import pickle
import cPickle
from numpy import *

verbose_utils = True
def claim_reset():
    if verbose_utils:
        print '''Resetting:
   '''+inspect.stack()[1][3]+'''
'''

#if writenet is called without a name, use the value
#default.
def writenet(name = 'default', value = None, hardcopy = True,np = False):
    caller_name = inspect.stack()[1][3]
    savename = caller_name + '_' +name
 
    globals()['lastname_'+caller_name] = name
    globals()['last_'+caller_name] = value

    if hardcopy:
        path = os.path.join(os.environ['NETWORK_TEMPPATH'], savename)
        f = open(path,'w')
        if np:
            save( f,value)
        else:
            pickle.dump(value, f)
        f.close()

def wn2(name = 'default',value = None, hardcopy = True, np = False):
    caller_name = inspect.stack()[1][3]
    savename = caller_name + '_' +name
 
    globals()['lastname_'+caller_name] = name
    globals()['last_'+caller_name] = value

    if hardcopy:
        path = os.path.join(os.environ['NETWORK_TEMPPATH'], savename)
        f = open(path,'w')
        if np:
            save( f,value)
        else:
            pickle.dump(value, f)
        f.close()


#if readnet is called without a name, it just
#returns the most recently globalized net.
def rn2(name = 'last', hardcopy = True,np = False):
    sxs = False
    out = -1

    caller_name = inspect.stack()[1][3]
    try:
        lname = globals()['lastname_'+caller_name]
        if lname == name or name == 'last':
            return globals()['last_'+caller_name], True
        else: 
            raise Exception()
    except Exception as e:
        if not hardcopy:
            return out, sxs
        savename = caller_name +'_'+name
        path =  os.path.join(os.environ['NETWORK_TEMPPATH'],savename)

        if hardcopy:
            path = os.path.join(os.environ['NETWORK_TEMPPATH'], savename)
            f = open(path,'r')
            if np:
                out = load(f)
                sxs = True
            else:
                out = pickle.load( f)
                sxs = True
            f.close()

    
        globals()['lastname_'+caller_name] = name
        globals()['last_'+caller_name] = out

    return out, sxs


#if readnet is called without a name, it just
#returns the most recently globalized net.
def readnet(name = 'last', hardcopy = True,np = False):
    sxs = False
    out = -1

    caller_name = inspect.stack()[1][3]
    try:
        lname = globals()['lastname_'+caller_name]
        if lname == name or name == 'last':
            return globals()['last_'+caller_name], True
        else: 
            raise Exception()
    except Exception as e:
        if not hardcopy:
            return out, sxs
        savename = caller_name +'_'+name
        path =  os.path.join(os.environ['NETWORK_TEMPPATH'],savename)

        if hardcopy:
            path = os.path.join(os.environ['NETWORK_TEMPPATH'], savename)
            f = open(path,'r')
            if np:
                out = load(f)
                sxs = True
            else:
                out = pickle.load( f)
                sxs = True
            f.close()

    
        globals()['lastname_'+caller_name] = name
        globals()['last_'+caller_name] = out
    if not sxs:
        raise Exception('readnet failed')
    return out


