#!/usr/bin/env python

import compbio.config as config
import compbio.broad.make_bsubs as make_bsubs
import sys, inspect, os, pickle
from Bio import AlignIO

def make_params():

  inp_dicts =[]
  for r, d, fs in os.walk(config.dataPath('::unseen_data')):
    for f in fs:
      print fs
      if '.stk' in f:
        filename = os.path.join(r,f)
        file_url = '::'+filename[filename.index('unseen'):]
        inp_dicts.append({'file_url':file_url})
        print f
      
  make_bsubs.make('process_stks.py', inp_dicts, mem_req = 400)

def main():

  name= sys.argv[1]

  f0 = open(os.environ['HOME'] + '/hello.txt','w')
  f0.write(name.__str__())
  f0.close()

  inp_dic = pickle.load(open(config.scriptInputPath(name)))
  
  filename = inp_dic['file_url']
  align = AlignIO.parse(open(config.dataPath(file_url)), 'stockholm')
  r0 = align.next()
  print 'hey whatsup'
  

  fout = open(config.scriptOutputPath(name), 'w')
  fout.write(r0.__str__())
  fout.close()

if __name__ == '__main__':
  main()
  exit()
