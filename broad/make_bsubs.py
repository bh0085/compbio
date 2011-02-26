import compbio.config as config
import subprocess
import os
import pickle

def make(script_name,input_dicts, mem_req = None):
  scr_path = os.path.join(config.root, 'scripts')
  bsub_path = os.path.join(scr_path, script_name + '.bsub')
  fopen = open(bsub_path, 'w')
    
  
  for i in range(len(input_dicts)):
    d = input_dicts[i]
    inp_file = os.path.join(scr_path, 'scr_inputs/'+ script_name+ str(i))
    inp_handle = open(inp_file,'w')
    out_file = os.path.join(scr_path, 'scr_output/'+ script_name+ str(i))
    pickle.dump(d,inp_handle)
    
    l = 'bsub -q compbio-week %s -i %s -o %s %s' % (\
      (lambda x: x == None and ' ' or ' -R %s ' % (x))(mem_req),
      inp_file, out_file, os.path.join(scr_path, script_name))
    # subprocess.Popen(l, shell = True)
    print l

    
