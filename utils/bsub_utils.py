#!/usr/bin/env python

import sys  
import subprocess, re
from numpy import roll
import simplejson as sjson

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
    return load_out(run_id)



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

def save_inp(inp_dict,run_id):
  input_dir = cfg.dataPath('batch/inputs')
  input_name = os.path.join(input_dir, run_id + '.inp')
  inp_file = open(input_name, 'w')
  pickle.dump(inp_dict, inp_file)
  inp_file.close()
  return input_name

def load_inp(run_id):
  input_dir = cfg.dataPath('batch/inputs')
  input_name = os.path.join(input_dir, run_id + '.inp')
  inp_file = open(input_name)
  inp_dict = pickle.load( inp_file)
  inp_file.close()
  return(inp_dict)


def save_out(out_dict, run_id):
  out_dir = cfg.dataPath('batch/outputs')
  out_name = os.path.join(out_dir, run_id + '.out')
  out_file = open(out_name, 'w')
  pickle.dump(out_dict, out_file)
  out_file.close()
  return out_name

def load_out(run_id):
  input_dir = cfg.dataPath('batch/outputs')
  input_name = os.path.join(input_dir, run_id + '.out')
  if not os.path.isfile(input_name):
    return None
  inp_file = open(input_name)
  inp_dict = pickle.load( inp_file)
  inp_file.close()
  return(inp_dict)

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
    if sys.argv[1] == 'bjobs':
        ids = sys.argv[2:]
        sys.stdout.write(sjson.dumps(bjobs(ids)))
        exit(0)
    elif sys.argv[1] == 'bout':
        run_id = sys.argv[2]
        sys.stdout.write(sjson.dumps(bout(run_id)))
        exit(0)                  
    else:
        raise Exception()

