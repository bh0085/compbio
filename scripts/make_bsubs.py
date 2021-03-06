import compbio.config as config
import subprocess

def make_bsubs(script_name,input_dicts, mem_req = None):
  scr_path = os.path.join(config.root, 'scripts')
  bsub_path = os.path.join(scr_path, script_name + '.bsub')
  fopen = open(bsub_path, 'w')
    
  
  for i in range(len(input_dicts)):
    d = input_dicts[i]
    inp_file = os.path.join(scr_path, 'scr_inputs/'+ script_name+ str(d))
    inp_handle = open(inp_file)
    out_file = os.path.join(scr_path, 'scr_output/'+ script_name+ str(d))
    pickle.dump(d,inp_handle)
    
    l = 'bsub -q compbio-week {3} -i {0} -o {1} {2}'.format(\
      inp_file, out_file, os.path.join(scr_path, script_name),\
        (lambda x: x == None and ' ' or ' -R {0} '.format(x))(mem_req))
    subprocess.popen(l, shell = True)
    print l

    
