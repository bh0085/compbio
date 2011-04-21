import pickle, os, itertools as it,re
import compbio.config as cfg
import inspect, subprocess as spc, time
import compbio.utils.remote_utils as rutils
#import compbio.utils.bsub_utils as butils
from compbio.utils.bsub_utils import *
from numpy import *
import simplejson as sjson

'''
Two classes and a disorganized bunch of routines for running bsub.

local_launcher: allows spawning of single bsub tasks from the local machine.
eyeball:        allows the spawing of many bsub tasks

Both eyeball and local_launcher can check up on the status of spawned threads.

'''

for p in ['batch','batch/inputs','batch/outputs','batch/logs','batch/tmp','batch/eye']:
  if not os.path.isdir(cfg.dataPath(p)):
    os.mkdir(cfg.dataPath(p))

class local_launcher(object):
  def __init__(self, scriptfile, scriptroot,
               func = None, run_id = None, host = 'tin'):
    '''
When calling scripts remotely I still have not worked out the most
clear possible syntax. Right now, it is up to the user to specify
both the remote path of the script and the remote path of the log
directory which she gets most easily via the function config.remotePath().

Anyway, having gotten those this script is pretty straightforward:

1) Call cmd = bsub.cmd().
2) Shellquote cmd and prepend: 'ssh tin'.
3) Run the final command via ssh.

Right now, all this version of local launch does is call the function 
remote_make_tests which runs a batch of clustering algorithms in matlab.
'''
    self.host = host
    self.func = func
    self.scriptfile = scriptfile
    print 'Fetching remote script path'
    remote_script =  cfg.remotePath(scriptfile,
                                       root = scriptroot)
    print 'Fetching remote log path'
    remote_logpath=  cfg.remotePath(cfg.dataPath('batch/logs'))

    self.run_id = get_run_id(1, 'bcl_ll')
    print 'Creating bsub commands'
    self.bscmd = cmd(remote_script, 
                func = func, 
                run_id = self.run_id,
                log_dir =remote_logpath)

  def launch(self):
    print 'Launching job for {0}.{1}'.format(self.scriptfile, self.func)
    comm = rutils.ssh_exec(self.bscmd, host = self.host)
    self.run_jobid = re.compile('Job <([\d]+)>').\
        search(comm).group(1)
    print '  submitted with jobID: {0}'.format(self.run_jobid)

  def status(self):
    scrpath= '${COMPBIO_PATH}/utils/bsub_utils.py'
    cmd = '{0} {1} {2}'.format(scrpath, 'bjobs', self.run_jobid)
    return sjson.loads(rutils.ssh_exec(cmd, host =self.host))

  def output(self):
    scrpath =  '${COMPBIO_PATH}/utils/bsub_utils.py'
    cmd = '{0} {1} {2}'.format(scrpath, 'bout', self.run_id)
    return sjson.loads(rutils.ssh_exec(cmd, host =self.host))

class eyeball(object):
  '''The class eyeball wraps all calls to bsub in my libraries. 

Eyeball has methods to check the status of bsub jobs underway as well as to recover program output once runs have terminated.

'''
  def __init__(self,
               scr_path, inp_dicts,
               func = None, 
               datapath = 'batch/eye/last.out',
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

    #Set up internally useful vars.
    self.datapath = datapath
    self.run_jobids = []

    #Use 'name' or 'scr' to set a prefix for runid generation
    if name == None:
      runid_prefix = os.path.splitext(\
        os.path.basename(scr_path))[0][-5:]
    else:
      runid_prefix = name[-5:]
    self.name = runid_prefix

    #Make a list of commands and corresponding run_ids.
    cmds = []
    self.run_names , self.run_ids = [], []
    for idx,  d in enumerate(inp_dicts):
      run_id=get_run_id(idx, prefix = runid_prefix)
      save_inp(d, run_id)
      self.run_names.append(run_id)
      cmds.append(cmd(scr_path,
                      func = func,
                      run_id = run_id))
    self.cmds = cmds
  def launch(self):
    for c in self.cmds:
      out = spc.Popen(c, stdout = spc.PIPE, \
                        shell = True).\
                        communicate()[0]
      self.run_jobids.append(re.compile('Job <([\d]+)>').\
                            search(out).group(1))

  def statii(self):
    '''
    Return the run statuses of programs launched under the control of
this eye. Uses bjobs.
'''
    jobs = bjobs(self.run_jobids) 
    statii = [j['STAT'].strip() for j in jobs.values()]
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
    #is what that run_jobids field is named after
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
  def package(self):
    outputs = self.outputs()
    pickle.dump(outputs,
                open(cfg.dataPath(self.datapath),'w'))
  def export(self, where = 'gliese'):
    cmdstr ='''
echo "connecting to remote server"; 
ssh gliese '
echo "copying files to remote server"; 
scp tin:{1} `python -c '\\''import compbio.config as config ; print config.dataPath("{2}")'\\''`; 
echo "copying successful"; 
exit' 
'''.format(where,
           cfg.dataPath(self.datapath), 
           self.datapath)

    print cmdstr    
    stdout = spc.Popen( cmdstr, 
                               shell = True,
                               stdout = spc.PIPE).communicate()
    print stdout
    return
  def awaitAndExport(self):
    print '''Entering a loop to await completion of all tasks for eye with name {0}'''.format(self.name)
    count = 0
    while 1:
      count += 1
      print '''Looping. [iter = {0}]'''.format(count)
      print 'statii:'
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
      time.sleep(10)

      
    self.package()
    self.export()





def cmd(scr_path, *args, **kwargs ):
  '''
Get a bsub command running the given sript

The templated way to call bsub.cmd() assumes that 'func' 
is specified as a kwarg and that run_id is specified or 
inferred from [name,num].

In this case, cmd() return a bsub command that with func
and run_id as the first and second arguments of scr_path.
In this case, scr_path can call func with argument run_id
which may then request an input dictionary from bsub on 
demand.

inputs:
  scr_path: path to an executable script, perhaps from
            config.getRemote().
  args: script args that will be inserted after func, 
        run_id if they are available.

keywords:
  func:         name of the subroutine of scr_path to call.
  [run_id | name]: either a run_id or a prefix for automatic
                generation of a run_id.
  project:      a project to mark the bsub job with.
  mem: [1GB]    any custom memory requests
  

'''
  func, name, run_id,project , Rmem, log_dir=\
      [kwargs.get('func', None),
       kwargs.get('name', ''),
       kwargs.get('run_id',''),
       kwargs.get('project','default'),
       kwargs.get('mem', '1'), 
       kwargs.get('log_dir',cfg.dataPath('batch/logs'))]

  #Set up a run_id if none is given.
  #By default increments to the max runid found.
  if not run_id:
    num = get_run_num()
    if name: prefix = name[:5]
    else: prefix = blnk
    run_id = get_run_id(num, prefix = prefix)

  if func != None:
    run_str = scr_path+' '+' '.join([func, run_id]+list(args)) 
  else:
    run_str = scr_path+ ' '  + ' '.join(args) 
  cmd = 'bsub -q compbio-week -J {3} -o {2} -P {0} -R \'rusage[mem={4}]\'  "{1}" '\
      .format(project, 
              run_str,
              os.path.join(log_dir,'{0}.log'.format(run_id) ),
              run_id,
              Rmem)

  
  return cmd

