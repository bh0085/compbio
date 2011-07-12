#!/usr/bin/env python
import cb.config as cfg
import os
import cb.utils.bsub_utils as bsu
import cb.utils.bsub as bsub
import subprocess

sdp = cfg.dataPath('soheil')
matfiles = [os.path.join(sdp,f) for f in os.listdir(sdp) if '.mat' in f]

for f in matfiles: 
    args = \
        f, \
        os.path.join(*( os.path.split(f)[:-1]+('mat_out',)+os.path.split(f)[-1:])),\
        os.path.join(*( os.path.split(f)[:-1]+('mat_all_out',)+os.path.split(f)[-1:]))

    script = 'test_mine'
    mat_cmd = '''"{2}('{0}', '{1}' ); exit"'''.format(*(args + (script,)))
    cstr = '''echo {0} | matlab -nojvm -nodisplay -nosplash '''.format(mat_cmd)
    bscmd = bsub.cmd(cstr, run_id = 'CLUSTER_'+os.path.split(f)[-1][:-4])
    subprocess.Popen(bscmd, shell = True)
    
                     
    
    break

exit
