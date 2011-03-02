#!/usr/bin/env python

import compbio.config as config
import compbio.broad.make_bsubs as make_bsubs
import sys, inspect, os, pickle
from Bio import AlignIO
import compbio.projects.cbdb.gb_alignment as gba
def make_params():

  inp_dicts =[]
  for r, d, fs in os.walk(config.dataPath('::unseen_data')):
    for f in fs:
      print fs
      if '.stk' in f:
        filename = os.path.join(r,f)
        file_url = '::'+filename[filename.index('unseen'):]
        inp_dicts.append({'file_url':file_url,
                          'ali_name':os.path.splitext(os.path.basename(f))[0]})
        print f
      
  make_bsubs.make('process_stks.py', inp_dicts, mem_req = 2000)

def main():

  name= sys.argv[1]


  print 'getting input dictionary'
  inp_dic = pickle.load(open(config.scriptInputPath(name)))
  
  file_url = inp_dic['file_url']
  print 'file_url: %s' % file_url
  print 'parsing alignment'

  p = config.dataPath(file_url)
  gba.fill_from_rfam_stk( p, reset = True)
  
if __name__ == '__main__':
  main()
  exit()
