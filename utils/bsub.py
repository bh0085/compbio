import pickle, os, itertools as it,re
import compbio.config as config
for p in ['batch','batch/inputs','batch/outputs','batch/logs']:
  if not os.path.isdir(config.dataPath(p)):
    os.mkdir(config.dataPath(p))

class eyeball(object):
  def __init__(self,
               scr_path, scriptargs, inp_dicts,
               name = None
           ):
    
    if name == None:
      runid_prefix = os.path.basename(scr_path) + '_'.join([''] + scriptargs)
    else:
      runid_prefix = runid_prefix

    do_bsub = True
    cmds = []
    self.run_ids = []
    for idx,  d in enumerate(batch_pdicts):
      run_id=bsub.get_run_id(idx, prefix = runid_prefix)
      bsub.save_inp(d, run_id)
      self.run_ids.append(run_ids)
      cmds.append(bsub.cmd(os.path.abspath(inspect.stack()[0][1]),\
                             ' '.join(scriptargs),\
                             run_id,\
                             do_bsub = do_bsub,\
                             run_id = run_id))
    if run:
      for c in cmds:
        out = subprocess.Popen(c, stdout = subprocess.PIPE, shell = True).\
            communicate()
        

  def statii():
    jobs = subprocess.Popen'bjobs '+ ' '.join(self.run_ids), shell = True, stdout = subprocess.PIPE).\
        communicate()[0]
    lines = '\n'.split(jobs)
    cols, lines = lines[0],lines[1:]
    
    
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
    cmd = 'bsub -q compbio-week -o {2} -P {0} -R \'rusage[mem=3]\'  "{1}" '\
        .format(project, 
                run_str,
                os.path.join(config.dataPath('batch/logs'),\
                               '{0}.log'.format(run_id) ))
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
