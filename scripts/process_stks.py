#!/usr/bin/env python


import sys, inspect, os, pickle
from Bio import AlignIO

os.chdir(os.path.dirname(inspect.stack()[0][1]))
inp, out = sys.argv[1:]
inp_dic = pickle.load(open(os.path.join('scr_inputs', inp)))
filename = inp_dic['filename']
align = AlignIO.parse(filename)
r0 = align.next()

fout = open(out, 'w')
fout.write(r0.__str__())
fout.close()

exit()
