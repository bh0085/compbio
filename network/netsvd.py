import pickle
from numpy import *
import matplotlib.pyplot as plt
import netutils as nu
import numpy as np
import numpy.linalg as lin
import netwriter as nw

def net_square_svd(name = nu.default_name, reset = 0):
    hardcopy = False
    try:
        if reset: raise Exception('compute')
        return nw.readnet(name, hardcopy = hardcopy)

    except Exception as e:
        if e.args[0] != 'compute': raise Exception()
        nw.claim_reset()
        sqa = nu.net_square_affinity( reset = mod(reset,2))

        U, S, Vh = lin.svd(sqa)
        V = Vh.T
        uvecs = U[:,:len(S)]
        vvecs = V
        dosave = (sqa[1],(uvecs, S, vvecs))
        nw.writenet(name, dosave,hardcopy = hardcopy)
        print '''
Ran net square svd for the current values.

Nothing has yet been saved to HD.

'''
        if reset:
            net_square_svd_U(reset = 2)
            net_square_svd_V(reset = 2)
            net_square_svd_S(reset = 2)
            
            print '''
...saved hardcopies of U,S,V
'''
        return dosave

def net_square_svd_U(name = nu.default_name,reset = 0):
    hardcopy = True
    try:
        if reset: raise Exception('compute')
        return nw.readnet(name, hardcopy = hardcopy)

    except Exception as e:
        if e.args[0] != 'compute': raise Exception()
        svd = net_square_svd(name, reset = mod(reset,2)) 
        dosave = svd[1][0]
        nw.writenet(name, dosave,hardcopy = hardcopy, np = True)
        return dosave    


def net_square_svd_V(name = nu.default_name,reset = 0):
    hardcopy = True
    try:
        if reset: raise Exception('compute')
        return nw.readnet(name, hardcopy = hardcopy,np = True)

    except Exception as e:
        if e.args[0] != 'compute': raise Exception()
        svd = net_square_svd(name, reset = mod(reset,2)) 
        dosave = svd[1][2]
        nw.writenet(name, dosave,hardcopy = hardcopy,np = True)
        return dosave    

def net_square_svd_S(name = nu.default_name,reset = 0):
    hardcopy = True
    try:
        if reset: raise Exception('compute')
        return nw.readnet(name, hardcopy = hardcopy,np = True)

    except Exception as e:
        if e.args[0] != 'compute': raise Exception()
        svd = net_square_svd(name, reset = mod(reset,2)) 
        dosave = svd[1][1]
        nw.writenet(name, dosave,hardcopy = hardcopy,np = True)
        return dosave    

#NOTE: no key is being saved with ggsvd...
#will have to recall net_square_affinity
#to get it.
def net_ggsvd(name = nu.default_name, reset = 0 ):
    hardcopy = False
    try:
        if reset: raise Exception('compute')
        return nw.readnet(name, hardcopy = hardcopy)

    except Exception as e:
        if e.args[0] != 'compute': raise Exception()
        nw.claim_reset()
        ggn = nu.net_genegene_norm( reset = mod(reset,2))

        
        U, S, Vh = lin.svd(ggn)
        V = Vh.T
        dosave =  (U,S,V)
        nw.writenet(name, dosave,hardcopy = hardcopy)

        print '''
Ran ggsvd.

No hardcopies have been saved yet'''
        if reset:
            net_ggsvd_U(reset = 2)
            net_ggsvd_V(reset = 2)
            net_ggsvd_S(reset = 2)
            print '''
...saved hardcopies of U,S,V
'''
        return dosave

def net_ggsvd_U(name = nu.default_name,reset = 0):
    hardcopy = True
    try:
        if reset: raise Exception('compute')
        return nw.readnet(name, hardcopy = hardcopy,np = True)

    except Exception as e:
        if e.args[0] != 'compute': raise Exception()
        svd = net_ggsvd(name, reset = mod(reset,2)) 
        dosave = svd[0]
        nw.writenet(name, dosave,hardcopy = hardcopy,np = True)
        return dosave    

def net_ggsvd_V(name = nu.default_name,reset = 0):
    hardcopy = True
    try:
        if reset: raise Exception('compute')
        return nw.readnet(name, hardcopy = hardcopy,np = True)

    except Exception as e:
        if e.args[0] != 'compute': raise Exception()
        svd = net_ggsvd(name, reset = mod(reset,2)) 
        dosave = svd[2]
        nw.writenet(name, dosave,hardcopy = hardcopy,np = True)
        return dosave  

def net_ggsvd_S(name = nu.default_name,reset = 0):
    hardcopy = True
    try:
        if reset: raise Exception('compute')
        return nw.readnet(name, hardcopy = hardcopy,np = True)

    except Exception as e:
        if e.args[0] != 'compute': raise Exception()
        svd = net_ggsvd(name, reset = mod(reset,2)) 
        dosave = svd[1]
        nw.writenet(name, dosave,hardcopy = hardcopy,np = True)
        return dosave  

def net_svd(name = nu.default_name,reset = 0):
    hardcopy = True
    try:
        if reset: raise Exception('compute')
        return nw.readnet(name = nu.default_name, hardcopy = hardcopy)
    except Exception as e:
        if e.args[0] != 'compute': raise Exception()
        nw.claim_reset()
        a = nu.net_affinity( reset = mod(reset,2))[0]
        N = a[0]
        U, S, Vh = lin.svd(N)
        V = Vh.T
        U = U[:,len(S)]
        dosave = (a[1],(U,S,V))
        nw.writenet(name,dosave,hardcopy = hardcopy)
        return dosave
