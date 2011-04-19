import pickle, os, itertools as it,re
import compbio.config as config
import inspect, subprocess

for p in ['batch','batch/inputs','batch/outputs','batch/logs']:
  if not os.path.isdir(config.dataPath(p)):
    os.mkdir(config.dataPath(p))

class eyeball(object):
  def __init__(self,
               scr_path, scriptargs, inp_dicts,
               name = None
           ):
    
    if name == None:
      runid_prefix = os.path.splitext(os.path.basename(scr_path))[0][-5:]
    else:
      runid_prefix = runid_prefix[-5:]


    do_bsub = True
    cmds = []
    self.run_names , self.run_ids = [], []
    for idx,  d in enumerate(inp_dicts):
      run_id=get_run_id(idx, prefix = runid_prefix)
      save_inp(d, run_id)
      self.run_names.append(run_id)
      cmds.append(cmd(scr_path,\
                        ' '.join(scriptargs),\
                        run_id,\
                        do_bsub = do_bsub,\
                        run_id = run_id))
    for c in cmds:
      out = subprocess.Popen(c, stdout = subprocess.PIPE, shell = True).\
          communicate()[0]
      self.run_ids.append(re.compile('Job <([\d]+)>').search(out).group(1))
    
        

  def statii(self):
    jobs = subprocess.Popen('bjobs '+ ' '.join(self.run_ids), shell = True, stdout = subprocess.PIPE).\
        communicate()[0]
    lines = jobs.split('\n')
    cols, lines = lines[0],lines[1:]
    for term in cols.split('\S+'):
      col_starts(term = cols.index(term))
    raise Exception()
    return lines
    
  def await():
    pass


def cmd(scr_path, *args, **kwargs ):

  run_id,project , inp_dict ,do_bsub =\
      [kwargs.get('run_id',''),
       kwargs.get('project','default'),
       kwargs.get('inp_dict',{}),
       kwargs.get('do_bsub', True)]
  assert run_id
  run_str = scr_path+ ' '  + ' '.join(args) 
  if do_bsub:
    cmd = 'bsub -q compbio-week -J {3} -o {2} -P {0} -R \'rusage[mem=3]\'  "{1}" '\
        .format(project, 
                run_str,
                os.path.join(config.dataPath('batch/logs'),\
                               '{0}.log'.format(run_id) ),
                run_id)
  else:
    cmd = '{0}'\
        .format(run_str)
  
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
  inp_file.close()

def load_inp(run_id):
  input_dir = config.dataPath('batch/inputs')
  input_name = os.path.join(input_dir, run_id + '.inp')
  inp_file = open(input_name)
  inp_dict = pickle.load( inp_file)
  inp_file.close()
  return(inp_dict)


def save_out(out_dict, run_id):
  out_dir = config.dataPath('batch/outputs')
  out_name = os.path.join(out_dir, run_id + '.out')
  out_file = open(out_name, 'w')
  pickle.dump(out_dict, out_file)
  out_file.close()
