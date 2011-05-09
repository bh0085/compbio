import scipy.io as sio
from numpy import *
import compbio.config as cfg
import matplotlib.pyplot as plt

def view():

    fnums = range(3000)[::5]
    #outputs =[ load_data( open(cfg.dataPath('batch/tmp/run_mcmc_{0:05}_tmp001.mat'.\
    #  
    outputs =[ sio.loadmat(cfg.dataPath('batch/tmp/mcmc_{0:05}_tmp001.mat'.format(num)))  for num in fnums ]

    douts = []
    for output in outputs:
        o00 = output['pout'][0][0]
        dout = dict([(k, o00[i]) for i , k in enumerate([elt[0] for elt in o00.dtype.descr])])
        douts.append(dout)


    raise Exception()
    plt.scatter( *array([(squeeze(o['stay_same']),squeeze(o['improve_ratio'])) for o in douts]).T )
    
    
    
