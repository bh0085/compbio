import pickle, os, itertools as it,re
import compbio.config as cfg
import inspect, subprocess as spc, time
import compbio.utils.remote_utils as rutils
from compbio.utils.bsub_utils import *
import compbio.utils.bsjobs
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
    self.scp_proc = None
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
                     log_dir =remote_logpath,
                     do_clear = True)

  def launch(self):
    print 'Launching job for {0}.{1}'.format(self.scriptfile, self.func)
    comm = rutils.ssh_exec(self.bscmd, host = self.host)
    self.run_jobid = re.compile('Job <([\d]+)>').\
        search(comm).group(1)
    print '  submitted with jobID: {0}'.format(self.run_jobid)

  def remote_status(self):
    scrpath= '${COMPBIO_PATH}/utils/bsruns.py'
    cmd = '{0} {1} {2}'.format(scrpath, 'bstatus', self.run_id)
    return sjson.loads(rutils.ssh_exec(cmd, host =self.host))
  def remote_output(self):
    '''NOT USED RIGHT NOW!
    This method relies on using stdout to communicate
    the output of remote execution but this seems silly.
    
    Calls that used to use this should migrate to remote_status
    which by convention deals only in variables short enough to 
    be transmitted by stdout
    '''

    scrpath =  '${COMPBIO_PATH}/utils/bsruns.py'
    cmd = '{0} {1} {2}'.format(scrpath, 'bout', self.run_id)
    return sjson.loads(rutils.ssh_exec(cmd, host =self.host))
  def fetch_start(self):
    '''
    Make sure this is not called until execution is finished.
    Otheriwse, the key 'outfile' may not yet be defined in the status dict.
    '''
    self.outfile = self.remote_status()['outfile']
    self.scp_proc = rutils.scp_data(self.outfile, self.outfile,
                               src_host = self.host)
  def fetch_await(self,maxtime = 10):
    print 'Checking status of the remote job manager.'
    print '...'
    while 1:
      status = self.remote_status()
      stat_str = status['status']
      print stat_str


      if stat_str in ['PEND' , 'RUN', 'UNK']:
        pass
      elif stat_str == 'EXIT':
        raise Exception()
      elif stat_str == 'DONE':
        break

      time.sleep(20)

      
      
    if not self.scp_proc: self.fetch_start()
    returnval = rutils.comm_timeout(self.scp_proc, 10)
    return returnval
  def output(self):
    '''
    Once the output has been fetched via scp from the remote
    process, load it up in the form of a python dictionary.
    '''
    fopen = open(cfg.dataPath(self.outfile))
    outval = pickle.load(fopen)
    return outval
    

  def quickRun(self):
    self.launch()
    self.fetch_await()
    return self.output()

class eyeball(object):
  '''The class eyeball wraps all calls to bsub in my libraries. 

Eyeball has methods to check the status of bsub jobs underway as well as to recover program output once runs have terminated.

'''
  def __init__(self,
               run_id,
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
    self.run_id = run_id

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
      save_data(d, run_id,'input')
      self.run_names.append(run_id)
      cmds.append(cmd(scr_path,
                      func = func,
                      run_id = run_id,
                      do_clear = True))
    self.cmds = cmds
  def launch(self):
    for idx, c in enumerate(self.cmds):
      out = spc.Popen(c, stdout = spc.PIPE, \
                        shell = True).\
                        communicate()[0]
      self.run_jobids.append(re.compile('Job <([\d]+)>').\
                            search(out).group(1))
      if mod(idx, 5) == 0 :
        save_data({'status':'RUN', 'str':'Jobs launched: {0}'.format(idx)}, self.run_id, 'status')
  def statii(self):
    '''
    Return the run statuses of programs launched under the control of
this eye. Uses bjobs.
'''
    jobs = compbio.utils.bsjobs.bjobs(self.run_jobids) 
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
    return [load_data(run_id, 'output') if statii[idx] == 'DONE' else None 
            for idx, run_id in enumerate(self.run_names)]


  def inputs(self):
    '''
    Get the inputs of programs run by this eye.

Returns the dictionary of inputs.
'''
    ids = self.run_names
    self.ins = []
    for i in ids:
      self.ins.append(load_data(i, 'input'))
    return self.ins
  
  def unfinished(self):
    stats = self.statii()
    return [s for s in stats if not s == 'DONE']
  def package(self):
    outputs = self.outputs()
    pickle.dump(outputs,
                open(cfg.dataPath(self.datapath),'w'))
  def export(self, where = 'gliese'):
    '''UNUSED This is the wrong way to do this.'''
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
  def await(self):
    count = 0
    while 1:
      count += 1
      stat_str =  '''Status: waiting\nLooping. [iter = {0}]\n'''.format(count)
      stat_str += 'statii:\n'

      jobs = compbio.utils.bsjobs.bjobs(self.run_jobids) 
      statii = [j['STAT'].strip() for j in jobs.values()]

      svals = dict(DONE= 1,
                   EXIT= -1,
                   RUN= 0,
                   PEND= 0)
      for k in svals.keys():
        stat_str +=  '   {1}:{0:02d}\n'.format(statii.count(k),k)
      vals = array([svals[k] for k in statii])
      save_data({'status': 'RUN','jobs':jobs}, self.run_id, 'status')
      if len(nonzero(not_equal(vals,1))[0]) == 0:
        break
      if len(nonzero(equal(vals,-1))[0]) > 0:
        save_data({'status':'EXIT','jobs':jobs}, self.run_id, 'status')
        raise Exception('Sorry but one of your scripts failed: {0}'.format(array(self.run_ids)[equal(vals,-1)]))
      time.sleep(20)
    save_data({'status':'RUN','jobs':jobs}, self.run_id, 'status')
                
  def complete(self):
    save_data({'status':'DONE', 'outfile': self.datapath},self.run_id, 'status')




def launch_actions(launch_instance, run_id):
  '''Mostly just clear out any current outputs and statuses.

  Unused right now.'''
  clear_data(run_id, ['output',  'status'])

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
  do_clear: [T/F] Clear any prior output/status msgs before running. 
  

'''
  func, name, run_id,project , Rmem, log_dir, do_clear=\
      [kwargs.get('func', None),
       kwargs.get('name', ''),
       kwargs.get('run_id',''),
       kwargs.get('project','default'),
       kwargs.get('mem', '1'), 
       kwargs.get('log_dir',cfg.dataPath('batch/logs')),
       kwargs.get('do_clear',True)]

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
  sub_cmd = 'bsub -q compbio-week -J {3} -o {2} -P {0} -R \'rusage[mem={4}]\'  "{1}" '\
      .format(project, 
              run_str,
              os.path.join(log_dir,'{0}.log'.format(run_id) ),
              run_id,
              Rmem)
  if do_clear:
    clr_cmd = '''${COMPBIO_PATH}'''+'''/utils/bsruns.py bclear {0}'''.\
        format(run_id)
    cmd = '; '.join([sub_cmd,clr_cmd])
  else:
    cmd = sub_cmd

  
  return cmd

