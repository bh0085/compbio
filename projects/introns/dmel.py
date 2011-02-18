from compbio import config
import os, re, itertools as it
from Bio import SeqIO

from compbio.biopython import biohilbert
from compbio.utils import memo as mem
from compbio.tests import drag_hilbert2 as dh2

import matplotlib.pyplot as plt
from numpy import *
import numpy as np

def run(reset= False, hash_r = False, dmel = False):
  if dmel:
    
    gpath = os.path.join(config.root, 'genomes/DMel')
  else:
    gpath = os.path.join(config.root,'genomes/Bacteria')
  flist = []
  glist = []
  
  for root, dirs, files in os.walk(gpath):
    for f in files:
      if re.search(re.compile('\.fna'), f):
        flist.append(os.path.join(root,f))
 
      if re.search(re.compile('\.gbk'), f):
        glist.append(os.path.join(root,f))


  fmax = 2
  flist = flist[:fmax]
  glist = glist[:fmax]
  if reset:
 
    fas = []
    gbs = []
    for f in flist:
      fas.append( list(SeqIO.parse(f, 'fasta'))[0])
    for g in glist:
      gbs.append( list(SeqIO.parse(g, 'genbank'))[0])
    mem.write('dmel', (fas,gbs))
  else:
    out, sxs = mem.read()
    fas, gbs = out

  if len(gbs) == len(fas): nc = len(gbs) 
  else: raise Exception()

  g0 = gbs[0]
  excount = 0 
  
  #types = map(lambda x: x.type, g0.features)
  types = sorted(g0.features, key  = lambda x: x.type)
  tcounts = {}
  for k, g in it.groupby(types, key = lambda x:x.type ):
    tcounts[k] = list(g)

  bh = biohilbert.BioHilbert(angle = 85,
                             lvls = 4)
  bh.initGB(g0,
            ndraw = 1000,
            score_only = True)

  if reset or hash_r :
    speedhash = bh.makeHash()
    mem.write('dmel', speedhash, register = 'speedhash')
  else:
    out, sxs = mem.read('dmel', register = 'speedhash')
    if sxs: speedhash = out
    else: raise Exception()
  bh.setHash(speedhash)
    
  actor_kwargs = {
    'skeleton':4,
    'skelwidth':1,
    'skelalpha':.8,
    'skeloff':0
    

  
  dh2.show_actor(bh, actor_kwargs = actor_kwargs) 

  
