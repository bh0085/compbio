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
               func = None, run_id = None, host = 'tin',
               input_dicts = []):
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

    #NOTE THAT FOR SOME REASON THIS ALL FAILS WHEN HOST 
    #IS SET TO NONE... THE SCRIPTS FAIL TO APPEND THE
    #DATAPATH
    print 'Fetching remote script path'
    remote_script =  cfg.remotePath(scriptfile,
                                    root = scriptroot,
                                    host = host)
    print 'Fetching remote log path'
    remote_logpath=  cfg.remotePath(cfg.dataPath('batch/logs'),
                                    host = host)

    self.run_id = run_id #get_run_id(1, 'bcl_ll')
    print 'Done!'
    print

    #now, transport the input_dicts
    print 'Saving input dictionaries on the remote host'
    input_path = locate_data(self.run_id, 'input')
    input_path = input_path[input_path.index('batch'):]
    save_data(input_dicts, self.run_id, 'input')
    print '...copying'
    comm = rutils.scp_data(input_path,input_path,dest_host = host).communicate()[0]
    print 'Done!'
    print

    log_file =  os.path.join(remote_logpath,'{0}.log'.format(self.run_id) )
    self.remote_logpath = log_file
    print 'Creating bsub commands'
    self.bscmd = cmd(remote_script, 
                     func = func, 
                     run_id = self.run_id,
                     log_file = log_file,
                     do_clear = True)
    
    print
    print 'Done!'
    print
    print 'Launcher is ready to launch in the background or await with "quickRun"'


      
    
  def launch(self):
    print 'Launching job for {0}.{1}'.format(self.scriptfile, self.func)
    comm = rutils.ssh_exec(self.bscmd, host = self.host)
    self.run_jobid = re.compile('Job <([\d]+)>').\
        search(comm).group(1)
    print '  submitted with jobID: {0}'.format(self.run_jobid)

    
  def fetch_logfile(self, child = None):
    import subprocess as spc

    import cb.utils.bsruns as bsruns
    if child != None:
      prc = spc.Popen('ssh {1} "bjobs -l {0}"'.format(child, self.host), stdout = spc.PIPE, shell = True); 
      stat =  prc.stdout.read();
      logpath = ''.join([l.strip() 
                         for l in re.compile('Output File[^<]*<([^>]*)',re.M + re.DOTALL)\
                           .search(stat).group(1).splitlines()])
                         
    else:
      logpath = self.remote_logpath
 
    prc = spc.Popen('ssh {1} "cat {0}"'.format(logpath, self.host), stdout = spc.PIPE, shell = True); 
    log =  prc.stdout.read();
    runs =[item for item in list(re.compile('^Sender: LSF.{10}', re.M).split(log)) if item.strip() ]
            
    results = [re.compile('Subject:(?P<subject>[^\n]*).*'+
                          'Started at (?P<start>[^\n]*).*'+
                          'Results reported at (?P<end>[^\n]*).*'+
                          'The output \(if any\) follows:(?P<output>.*)', re.M + re.DOTALL)
               .search(r).groupdict()
               for r in runs]
    return results
    
    

  def remote_status(self):
    scrpath= '${COMPBIO_PATH}/utils/bsruns.py'
    cmd = '{0} {1} {2}'.format(scrpath, 'bstatus', self.run_id)
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
               scr_path, 
               inp_dicts,
               mem = 2,
               func = None, 
               datapath = 'batch/eye/',
               name = None
           ):
    '''
Create an eye with a set of jobs. Run them.

NOTE 1: This was originally designed to construct
run_ids on the fly but now, if input dictionaries have
the field run_id defined, it will use the run_id value
as a run_id for each process.

inputs:
  scr_path:    the path of the script to be run. [String]
  scriptargs:  command line args to the script.  [String x nargs]
  inp_dicts:   input dicts for each instance of script
               implicitly sets the number of calls to script.

kwds:
  mem:         manually specify bsub memory requirements
  datapath:    if we choose to package the output after launch, 
               datapath specifies the location to save packaged data.
  name:        if run_ids are not specified, name is used to put a
               prefix on each run_id
  func:        func to run

'''

    #Set up internally useful vars.
    self.run_jobids = []
    self.run_id = run_id
    self.datapath = datapath + self.run_id + '.out'
    self.update_status('RUN', {'state':'beginning config'})

    #Use 'name' or 'scr' to set a prefix for runid generation
    if name == None:
      runid_prefix = os.path.splitext(\
        os.path.basename(scr_path))[0][-5:]
    else:
      runid_prefix = name
    self.name = runid_prefix

    #Make a list of commands and corresponding run_ids.
    cmds = []
    self.run_names , self.run_ids = [], []
    for idx,  d in enumerate(inp_dicts):
      if d.has_key('run_id'):
        run_id = d['run_id']
      else:
        run_id=get_run_id(idx, prefix = runid_prefix)

      save_data(d, run_id,'input')
      self.run_names.append(run_id)
      cmds.append(cmd(scr_path,
                      func = func,
                      run_id = run_id,
                      do_clear = True, 
                      mem = mem))
    resets = zeros(len(self.run_names))
    self.cmds = cmds
    self.update_status('RUN',{'state':'finished config; unlaunched'})

  def launch(self):
    self.update_status('RUN',{'state':'no jobs launched'})
    prc_q = []
    #SUBMIT THE JOBS IN BATCHES OF SIZE sub_parallell_count
    sub_parallell_count = 10
    for c in self.cmds:      
      prc_q.append(spc.Popen(c, stdout = spc.PIPE, shell = True))
      if len(prc_q) >= sub_parallell_count:
        self.run_jobids.extend([re.compile('Job <([\d]+)>').\
                                  search(p.stdout.read())\
                                  .group(1) for p in prc_q])
        prc_q = []
        self.update_status('RUN',{'state':'launching jobs, {0} launched'.format(len(self.run_jobids))})
        
    for idx in range(len(self.cmds)): 
      print 'job{0}:'.format(idx)
      print self.cmds[idx]

  def statii(self):
    '''
    Return the run statuses of programs launched under the control of
this eye. Uses bjobs.
'''
    jobs = compbio.utils.bsjobs.bjobs(self.run_jobids) 
    statii = dict([( k, j['STAT'].strip()) for k, j in jobs.iteritems()])
    return statii


  def update_status(self,stat_str, other = {}):
    '''Note that during the course of execution, the status of the
task manager in bjobs will necessarily be 'RUN'. Meanwhile,
the children can have various statii.

Thus, the bsub statii of the children are used because they
are informative (and child statii will necessarily be of a valid
bsub status type) whereas the eyeball's status is arbitrary.

'''
    status = {'status':stat_str}
    status['children'] = self.statii()
    status.update(other)
    save_data(status, self.run_id, 'status')

  def outputs(self):
    '''
    Get the outputs of programs run by this eye.

For programs that have been run, returns the output dictionary.
For programs that have not yet been completed, or that have failed,
returns a dictionary with the run id and the cause of failure.
'''
    outputs = []
    statii = self.statii()
    for  run_id in self.run_names:
      try:
        if statii[run_id] == 'DONE':
          data = load_data(run_id, 'output')
        else:
          data = dict(failure = 'job_stat({1}): {0}'.format(statii[run_id], run_id),
                      job = run_id)
      except Exception(), e:
        data = dict(failure = e,text = 'Exception in eyeball collection\n probably, outfile does not exist) step', job =run_id)
      outputs.append(data)
    return outputs


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
      if len(nonzero(less(vals,1))) == 0:
        self.update_status('RUN',{'awaiting':'finished'})
        break
      self.update_status('RUN', {'awaiting':'pending:{0}, running:{1}, done:{2}, failed:{3}'.\
                                   format(statii.count('PEND'),statii.count('RUN'), 
                                          statii.count('DONE'),statii.count('EXIT'))})
      time.sleep(20)

                
  def complete(self):
    self.update_status('DONE', {'outfile': self.datapath})




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
  
  print
  print
  func, name, run_id,project , Rmem, log_file, log_dir, do_clear=\
      [kwargs.get('func', None),
       kwargs.get('name', ''),
       kwargs.get('run_id',''),
       kwargs.get('project','default'),
       kwargs.get('mem', '1'), 
       kwargs.get('log_file', None),
       kwargs.get('log_dir',cfg.dataPath('batch/logs')),
       kwargs.get('do_clear',True)]

  #Set up a run_id if none is given.
  #By default increments to the max runid found.
  if not run_id:
    num = get_run_num()
    if name: prefix = name[:5]
    else: prefix = blnk
    run_id = get_run_id(num, prefix = prefix)


  if log_file == None:
    print 'WARNING:     Writing a log file path for some reason.'
    log_file = os.path.join(log_dir,'{0}.log'.format(run_id) )
  print '{1:40}{0}'.format(log_file, 'logfile:')

  if func != None:
    run_str = scr_path+' '+' '.join([func, run_id]+list(args)) 
  else:
    print 'WARNING:    No function specified. Using alternate form of bscmd where run_str = path + join(args)'
    run_str = scr_path+ ' '  + ' '.join(args)     
  print '{1:40}{0}'.format(run_str,'run_path:')

  sub_cmd = 'bsub -q compbio-week -J {3} -o {2} -P {0} -R \'rusage[mem={4}]\'  "{1}" '\
      .format(project, 
              run_str,
              log_file,
              run_id,
              Rmem)

  print '{1:40}{0}'.format(do_clear,'clear_prev:')
  if do_clear:
    clr_cmd = '''${COMPBIO_PATH}'''+'''/utils/bsruns.py bclear {0}'''.\
        format(run_id)
    cmd = '; '.join([clr_cmd,sub_cmd])
  else:
    cmd = sub_cmd

  print '{1:40}{0}'.format(cmd,'command: ')
  
  return cmd

