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
    
    l = 'bsub -q compbio-week %s -i %s -o %s %s' % (\
      (lambda x: x == None and ' ' or ' -R %s ' % (x))(mem_req),
      inp_file, '${HOME}/bsub_logs/'+inp_file, script_name)
    ls.append(l)
  cmd = "ssh tin '" + '; '.join(ls) + "'"
  subprocess.Popen(cmd, shell = True)
