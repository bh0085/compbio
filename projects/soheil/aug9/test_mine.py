#!/usr/bin/env python
import cb.config as cfg
import os
import cb.utils.bsub_utils as bsu
import cb.utils.bsub as bsub
import subprocess
import pipes 


sdp = cfg.dataPath('soheil')
matfiles = [os.path.join(sdp,f) for f in os.listdir(sdp) if '.mat' in f]

for control in [0,2,3]:
  for f in matfiles: 
      args = \
          f, \
          os.path.join(*( os.path.split(f)[:-1]+('mat_out_{0}'.format(control),)+os.path.split(f)[-1:])),\
          os.path.join(*( os.path.split(f)[:-1]+('mat_all_out_{0}'.format(control),)+os.path.split(f)[-1:]))
  
      for a in args: 
          if not os.path.isdir(os.path.dirname(a)):
              os.makedirs(os.path.dirname(a))
              print 'made directory: {0}'.format(os.path.dirname(a))
  
      
      script = 'test_mine'
      mat_cmd ='''\\"{3}('{0}', '{1}', '{2}', '{3}' ); exit\\"'''.format(\
          *(args + (script,)+(control,)))

      print mat_cmd
      cstr = '''echo {0} | matlab -nojvm -nodisplay -nosplash '''.\
                             format(mat_cmd)
      bscmd = bsub.cmd(cstr, run_id = 'CLUSTER_'+os.path.split(f)[-1][:-4], mem = 4)
  
      prc = subprocess.Popen(bscmd, shell = True)
      prc.communicate()
      #print bscmd

exit
