import compbio.config as config
import subprocess
import os
import pickle

def make(script_name,input_dicts, mem_req = None):
  scr_path = os.path.join(config.root, 'scripts')    
  
  ls = []
  for i in range(len(input_dicts)):
    d = input_dicts[i]
    inp_file =  script_name+ str(i)
    inp_handle = open(config.scriptInputPath(inp_file),'w')
    out_file =script_name+ str(i)
    pickle.dump(d,inp_handle)
    
    l = 'bsub -q compbio-week %s -o %s %s %s' % (\
      (lambda x: x == None and ' ' or ' -R %s ' % \
         (x))(mem_req), '${HOME}/bsub_logs/'+inp_file, script_name, inp_file)
    ls.append(l)
  cmd = "ssh tin '" + '; '.join(ls) + "'"

  prc0 = subprocess.Popen("rsync -rv ${COMPBIO_PATH}/scripts/scr_inputs/ tin:'${COMPBIO_PATH}/scripts/scr_inputs/'", shell = True, stdout = subprocess.PIPE)
  c0 = prc0.communicate()
  print 'RSYNCed inputs directory'
  prc1 = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE)
  print 'Sent batch commands out'
