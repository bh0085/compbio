#!/usr/bin/env python

import compbio.config as config
import compbio.broad.make_bsubs as make_bsubs
import sys, inspect, os, pickle
from Bio import AlignIO

def make_params():

  inp_dicts =[]
  for r, d, fs in os.walk(config.dataPath(config.dataURL('unseen_data'))):
    for f in fs:
      if '.stk' in f:
        filename = os.path.join(r, f)
        inp_dicts.append({'filename':filename})
      
  make_bsubs.make('process_stks', inp_dicts)

def main():

  os.chdir(os.path.dirname(inspect.stack()[0][1]))
  inp, out = sys.argv[1:]
  print inp, out
  inp_dic = pickle.load(open(os.path.join('scr_inputs', inp)))
  filename = inp_dic['filename']
  align = AlignIO.parse(open(filename), 'stockholm')
  r0 = align.next()

  fout = open(out, 'w')
  fout.write(r0.__str__())
  fout.close()

if __name__ == '__main__':
  main()
  exit()
