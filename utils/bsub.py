import pickle, os, itertools as it,re
import compbio.config as config
for p in ['batch','batch/inputs','batch/outputs','batch/logs']:
  if not os.path.isdir(config.dataPath(p)):
    os.mkdir(config.dataPath(p))

def cmd(scr_path, run_id,
        project = 'default', 
        inp_dict = {}, 
        do_bsub = True):
    
  if do_bsub:
    cmd = 'bsub -q compbio-week -o {3} -P {0} "{1} {2}" '\
        .format(project, 
                scr_path, 
                run_id,
                os.path.join(config.dataPath('batch/logs'),\
                               '{0}.log'.format(run_id) ))
  else:
    cmd = '{0} {1}'\
        .format(scr_path, run_id)
  
  return cmd

def get_run_num():
  cur_id = max([int(e) for e in re.findall(\
          re.compile('([0-9]+)'),' '.join(it.chain(*
              [os.listdir(config.dataPath(d)) 
               for d in  ['batch/inputs',
                          'batch/outputs',
                          'batch/logs']])))]+ [-1])
  num = cur_id + 1
  return num
def get_run_id(num, prefix = 'R' ):
  return '{0}%05i'.format(prefix) % (num,)

def save_inp(inp_dict,run_id):
  input_dir = config.dataPath('batch/inputs')
  input_name = os.path.join(input_dir, run_id + '.inp')
  inp_file = open(input_name, 'w')
  pickle.dump(inp_dict, inp_file)
  print input_name
  inp_file.close()

def load_inp(run_id):
  input_dir = config.dataPath('batch/inputs')
  input_name = os.path.join(input_dir, run_id + '.inp')
  inp_file = open(input_name)
  inp_dict = pickle.load( inp_file)
  inp_file.close()
  return(inp_dict)
