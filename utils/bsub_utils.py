#!/usr/bin/env python
'''A bunch of utilities designed to get information about bsub jobs 
launched via the classes in bsub.py.

Functions are designed to be called within python. 
They mostly take run_ids as input. 
They return python objects

For command line executable scrips, look in 
  bsruns.py     (Scripts called taking run_ids as input)
  bsjobs.py     (Scripts called taking jobIDs as input)

Instead of returning python objects, scripts executed from the
commandline print json formatted strings to stdout. They can be
used for example to communicate data over ssh calls.
'''


import subprocess, re, os, sys, pickle
from numpy import roll
import simplejson as sjson
import compbio.config as cfg

#Utility functions
def get_run_num():
  '''
  Automatically get a run number from the max of all files
  so far saved in input/output/logs
  '''
  cur_id = max([int(e) for e in re.findall(\
          re.compile('([0-9]+)'),' '.join(it.chain(*
              [os.listdir(cfg.dataPath(d)) 
               for d in  ['batch/inputs',
                          'batch/outputs',
                          'batch/logs']])))]+ [-1])
  num = cur_id + 1
  return num
def get_run_id(num, prefix = 'R' ):
    '''
    Get a generic run_id.
    '''
  return '{0}%05i'.format(prefix) % (num,)


def locate_data(run_id, data_type):
    '''
    Find the savepath for a run to output/input a given datatype.
    Datatypes so far implemented are 'output', 'input', 'status'.
    
    'output' and 'input' types are stored in python serialized pickle
    format and may be arbitrarily large.

    'status' on the the other hand is stored in json formatted strings
    so that they can be read and communicated over ssh quickly, simply
    and without worry about buffering/ illegal character issues.

    usage: locate_data(run_id [String] , data_type [String] )

    '''

    templates = {'output':'batch/outputs/{0}.pickle',
                 'input':'batch/inputs/{0}.pickle',
                 'status':'batch/statii/{0}.json'}
    return cfg.dataPath(templates[data_type].format(run_id))
def clear_data(run_id, data_types):
    '''
    Remove run_id's data from the disc.

    input:
      run_id:     run_id with which to locate data to be removed
      data_types: a list of datatypes to be removed.
    '''
    for d in data_types:
        path = locate_data(run_id, data_type)
        if os.path.isfile(path):
            os.remove(path)
def save_data(contents, run_id, data_type):
    '''
    Save data at a location determined by run_id via
    the function locate data which looks up a path for
    a given datatype in its templates dict.

    datatypes include: 'input', 'output', 'status'
    '''
    path = locate_data(run_id, data_type)
    fopen = open(path, 'w')
    if path[-4:] == 'ckle':
        pickle.dump(contents, fopen)
    else:
        sjson.dump(contents, fopen)
    fopen.close()
def load_data(run_id, data_type):
    '''
    Load data at a location determined by run_id via
    the function locate data which looks up a path for
    a given datatype in its templates dict.

    datatypes include: 'input', 'output', 'status'
    '''
    path = locate_data(run_id, data_type)
    fopen = open(path)
    if path[-4:] == 'ckle':
        contents = pickle.load( fopen)
    else:
        contents = sjson.load( fopen)
    fopen.close()    
    return contents


def tmp_fnames(run_id, num):
  '''get temporary filenames that scripts can write to'''
  tmp_dir = cfg.dataPath('batch/tmp')
  names = [os.path.join(tmp_dir, run_id + '_tmp{0:03d}'.format(idx))
                        for idx in range(num)]
                   
  return names

def mat_tmp_fnames(run_id, num):
  '''get temporary filenames with .mat appended that matlab
  saves can be written to. 
  
  (matlab doesn't like loading save files without .mat extension).
  '''
  tmp_dir = cfg.dataPath('batch/tmp')
  names = [os.path.join(tmp_dir, run_id + '_tmp{0:03d}.mat'.format(idx))
                        for idx in range(num)]
                   
  return names
