import subprocess as spc
import compbio.utils.bsub as bsub
import scipy.io as sio
from numpy import *

def runmat(script, input_dict, run_id):
    tmpnames = bsub.mat_tmp_fnames(run_id,2)
    sio.savemat(tmpnames[0],input_dict)
    mat_cmd = '''"{2}('{0}', '{1}' ); exit"'''.format(*(tmpnames + [script]))
    cstr = '''echo {0} | matlab -nojvm -nodisplay -nosplash '''.format(mat_cmd)
    sxs = spc.Popen(cstr, shell = True).communicate()
    output = sio.loadmat(tmpnames[1])
    o00 = output['out_struct'][0][0]
    dout   = dict([(k, o00.__getattribute__(k)) for k in o00._fieldnames])
    return dout

    
