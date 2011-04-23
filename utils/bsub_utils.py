#!/usr/bin/env python

import subprocess, re, os, sys, pickle
from numpy import roll
import simplejson as sjson
import compbio.config as cfg

#Scripts designed to be called from the command line.
def bjobs(jobids):
    '''
    Return the run statuses of programs launched under the control of
this eye. Uses bjobs.
'''
    #Get active jobs.
    jobs = subprocess.Popen('bjobs '+ ' '.join(jobids), 
                            shell = True, 
                            stdout = subprocess.PIPE).\
                            communicate()[0]
    lines = jobs.split('\n')
    cols, lines = lines[0],lines[1:]
    col_starts = {}

    #Parse job list into columns
    terms =  re.compile('\s+').split(cols)
    starts = [cols.index(t) for t in terms]
    ends =   [cols.index(t) for t in roll(terms, -1)]
    ends[-1] = len(cols)
    col_ranges = dict([(t, (starts[i], ends[i])) 
                       for i, t in enumerate(terms)])

    job_dicts = {}
    for l in lines:
      if l.strip() == '': continue
      d0 = dict([(k, l[col_ranges[k][0]:col_ranges[k][1]]) 
                 for k in terms])
      run_id = d0['JOBID']
      job_dicts[run_id] = d0
    return(job_dicts)
def bout(run_id):
    return load_data(run_id, 'output')

def bstatus(run_id):
    '''Get the status of a currently running job.
By convention, statuses will be saved in json
format so that bstatus is compatible with printing
to the stdout and communication over ssh.

usage:       bstatus(run_id)
commandline: bsub_utils.py bstatus run_id
'''
    return load_data(run_id, 'status')

#Utility functions
def get_run_num():
  cur_id = max([int(e) for e in re.findall(\
          re.compile('([0-9]+)'),' '.join(it.chain(*
              [os.listdir(cfg.dataPath(d)) 
               for d in  ['batch/inputs',
                          'batch/outputs',
                          'batch/logs']])))]+ [-1])
  num = cur_id + 1
  return num
def get_run_id(num, prefix = 'R' ):
  return '{0}%05i'.format(prefix) % (num,)




def locate_data(run_id, data_type):
    templates = {'output':'batch/outputs/{0}.pickle',
                 'input':'batch/inputs/{0}.pickle',
                 'status':'batch/statii/{0}.json'}
    return cfg.dataPath(templates[data_type].format(run_id))
def clear_data(run_id, data_type):
    for d in data_types:
        path = locate_data(run_id, data_type)
        if os.path.isfile(path):
            os.remove(path)
def save_data(contents, run_id, data_type):
    path = locate_data(run_id, data_type)
    fopen = open(path, 'w')
    if path[-4:] == 'ckle':
        pickle.dump(contents, fopen)
    else:
        sjson.dump(contents, fopen)
    fopen.close()
def load_data(run_id, data_type):
    path = locate_data(run_id, data_type)
    fopen = open(path)
    if path[-4:] == 'ckle':
        contents = pickle.load( fopen)
    else:
        contents = sjson.load( fopen)
    fopen.close()    
    return contents
    
    



def tmp_fnames(run_id, num):
  tmp_dir = cfg.dataPath('batch/tmp')
  names = [os.path.join(tmp_dir, run_id + '_tmp{0:03d}'.format(idx))
                        for idx in range(num)]
                   
  return names

def mat_tmp_fnames(run_id, num):
  tmp_dir = cfg.dataPath('batch/tmp')
  names = [os.path.join(tmp_dir, run_id + '_tmp{0:03d}.mat'.format(idx))
                        for idx in range(num)]
                   
  return names

if __name__ == '__main__':
    assert len(sys.argv) > 2
    if sys.argv[1] in [ 'bjobs', 'bout', 'bstatus' ]:
        #Dump data to json and write to stdout
        ids = sys.argv[2:]
        sys.stdout.write(sjson.dumps(globals()[sys.argv[1]](ids)))
        exit(0)
    else:
        raise Exception()

