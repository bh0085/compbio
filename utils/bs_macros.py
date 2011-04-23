import subprocess as spc
import compbio.utils.bsub as bsub
import scipy.io as sio
from numpy import *

def runmat(script, input_dict, run_id):
    tmpnames = bsub.mat_tmp_fnames(run_id,2)
    sio.savemat(tmpnames[0],input_dict)

    cstr = '''matlab -nodesktop -nosplash -r "{2}('{0}', '{1}' )"'''
