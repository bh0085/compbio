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
def writenet(name = 'default', value = None, hardcopy = False,np = False):
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
def readnet(name = 'last', hardcopy = False,np = False):
    caller_name = inspect.stack()[1][3]
    print 'reading'
    try:
        lname = globals()['lastname_'+caller_name]
        print lname
        
        if lname == name or name == 'last':
            return globals()['last_'+caller_name]
        else: 
            raise Exception()
    except Exception as e:
        print e
        if not hardcopy:
            raise Exception('compute')
        savename = caller_name +'_'+name
        path =  os.path.join(os.environ['NETWORK_TEMPPATH'],savename)
        try:
                if hardcopy:
                    path = os.path.join(os.environ['NETWORK_TEMPPATH'], savename)
                    f = open(path,'w')
                    if np:
                        load(f)
                    else:
                        pickle.load( f)
                    f.close()
        
        except Exception as e2:
            print e2
            raise Exception('readnet failed: '+ name + ' unsaved for ' + caller_name)
        globals()['lastname_'+caller_name] = name
        globals()['last_'+caller_name] = value
        return value


