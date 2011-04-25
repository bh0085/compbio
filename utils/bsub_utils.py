#!/usr/bin/env python

import subprocess, re, os, sys, pickle
from numpy import roll
import simplejson as sjson
import compbio.config as cfg

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
