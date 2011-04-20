import subprocess as spc
import compbio.utils.bsub as bsub
import scipy.io as sio
from numpy import *

def runmat(script, input_dict, run_id):
    tmpnames = bsub.mat_tmp_fnames(run_id,2)
    sio.savemat(tmpnames[0],input_dict)
    cstr = '''matlab -nodesktop -nosplash -r "{2}('{0}', '{1}'''.\
        format(*(tmpnames+[script]) )
    sub = spc.Popen(cstr,shell = True, stdout = spc.PIPE).\
        communicate()[0]

    
    struct = sio.loadmat(tmpnames[1])['out_struct']
    o00= struct[0][0]
    fns = o00._fieldnames
    out_dict = dict( [(fn, o00.__getattribute__(fn)) 
                     for fn in fns ])
    return out_dict
