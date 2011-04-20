import pickle, os, itertools as it,re
import compbio.config as config
import inspect, subprocess, time

#A class and a bunch of routines for 
#running and tracking the results of bsub.
#
#The main visible class here is eyeball which
#can be called bit a script path and args as well
#as a list of dictionaries representing input.

for p in ['batch','batch/inputs','batch/outputs','batch/logs','batch/tmp','batch/eye']:
  if not os.path.isdir(config.dataPath(p)):
    os.mkdir(config.dataPath(p))

class eyeball(object):
  '''The class eyeball wraps all calls to bsub in my libraries. 

Eyeball has methods to check the status of bsub jobs underway as well as to recover program output once runs have terminated.

'''
  def __init__(self,
               scr_path, scriptargs, inp_dicts,
               data = 'batch/eye/last.out',
               name = None
           ):
    '''
    Create an eye with a set of jobs. Run them.

inputs:
  scr_path:    the path of the script to be run. [String]
  scriptargs:  command line args to the script.  [String x nargs]
  inp_dicts:   input dicts for each instance of script
               implicitly sets the number of calls to script.

'''
    self.datapath = datapath
    if name == None:
      runid_prefix = os.path.splitext(os.path.basename(scr_path))[0][-5:]
    else:
      runid_prefix = runid_prefix[-5:]

    self.name = runid_prefix

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
    '''
    Return the run statuses of programs launched under the control of
this eye. Uses bjobs.
'''
    jobs = subprocess.Popen('bjobs '+ ' '.join(self.run_ids), shell = True, stdout = subprocess.PIPE).\
        communicate()[0]
    lines = jobs.split('\n')
    cols, lines = lines[0],lines[1:]
    col_starts = {}
    for term in re.compile('\s+').split(cols):
      col_starts.update([(term, cols.index(term))])
    
    statii = []
    for l in lines:
      statii.append( l[col_starts['STAT']:col_starts['QUEUE']].strip())
    statii = [s for s in statii if s != '']
    return statii
      
  def outputs(self):
    '''
    Get the outputs of programs run by this eye.

For programs that have been run, returns the output dictionary.
For programs that have not yet been completed, returns: None
'''
    #Note that even though the functions I am calling ask for 'ids',
    #I use the field i.run_names...
    #this is because I use a different id that than the bsub job id which
    #is what that run_ids field is named after
    statii = self.statii()
    return [load_out(run_id) if statii[idx] == 'DONE' else None 
            for idx, run_id in enumerate(self.run_names)]


  def inputs(self):
    '''
    Get the inputs of programs run by this eye.

Returns the dictionary of inputs.
'''
    ids = self.run_names
    self.ins = []
    for i in ids:
      self.ins.append(load_inp(i))
    return self.ins
  
  def unfinished(self):
    stats = self.statii()
    return [s for s in stats if not s == 'DONE']
  def package(self, out_datapath):
    outputs = self.outputs()
    pickle.save(open(config.dataPath(self.datapath),'w'),
                outputs)
  def export(self, where = 'gliese'):
    cmdstr ='''echo "connecting to remote server"; ssh gliese 'echo "copying files to remote server"; scp tin:{1} `python -c '\''import compbio.config as config ; print config.dataPath("{2}")'\''`; echo "copying successful"; exit' '''.format(where,config.dataPath(self.datapath), self.datapath)

    stdout = subprocess.Popen( cmdstr, shell = True,
                               stdout = subprocess.PIPE).communicate()
    print stdout
    return
  def awaitAndExport(self):
    print '''Entering a loop to await completion of all tasks for eye with name {0}'''.format(self.name)
    count = 0 
    while 1:
      count += 1
      print '''Looping. [iter = {0}]'''.format(count)
      print 'statii:'
      time.sleep(1)
      statii = self.statii()
      svals = dict(DONE= 1,
                   EXIT= -1,
                   RUN= 0,
                   PEND= 0)
      for k in svals.keys():
        print '   {1}:{0:02d}'.format(statii.count(k),k)
      vals = array([svals[k] for k in statii])
      if len(nonzero(not_equal(vals,1))[0]) == 0:
        break
      if len(nonzero(equal(vals,-1))[0]) > 0:
        raise Exception('Sorry but one of your scripts failed.')

    self.package()
    self.export()

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
  return input_name

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
  return out_name

def load_out(run_id):
  input_dir = config.dataPath('batch/outputs')
  input_name = os.path.join(input_dir, run_id + '.out')
  inp_file = open(input_name)
  inp_dict = pickle.load( inp_file)
  inp_file.close()
  return(inp_dict)

def tmp_fnames(run_id, num):
  tmp_dir = config.dataPath('batch/tmp')
  names = [os.path.join(tmp_dir, run_id + '_tmp{0:03d}'.format(idx))
                        for idx in range(num)]
                   
  return names

def mat_tmp_fnames(run_id, num):
  tmp_dir = config.dataPath('batch/tmp')
  names = [os.path.join(tmp_dir, run_id + '_tmp{0:03d}.mat'.format(idx))
                        for idx in range(num)]
                   
  return names
