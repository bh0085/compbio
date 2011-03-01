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
        filename = os.path.join(r, f)
        inp_dicts.append({'filename':filename})
        print f
      
  make_bsubs.make('process_stks', inp_dicts)

def main():

  inp, out = sys.argv[1:]

  inp_dic = pickle.load(open(config.scriptInputPath(inp)))
  filename = inp_dic['filename']
  align = AlignIO.parse(open(filename), 'stockholm')
  r0 = align.next()

  fout = open(config.scriptOutputPath(out), 'w')
  fout.write(r0.__str__())
  fout.close()

if __name__ == '__main__':
  main()
  exit()
