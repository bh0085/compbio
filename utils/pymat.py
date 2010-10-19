import os
import scipy.io as sio
import subprocess

def runmat(command, inp = {}):
    path = os.environ['PYMAT_PATH']
    fname = os.path.join(path, 'python.mat')
    f = (open(fname,'w'))
    sio.savemat(f,inp)
    f.close()

    ex_str = command + ''
    call_list = ['matlab','-nosplash','-nodesktop','-r '+ex_str]
                
    subprocess.call(call_list)
    outname = os.path.join(path,'matlab.mat')

    f = open( os.path.join(path , 'matlab.mat'))
    output = sio.loadmat(f)
    f.close()
    return output
