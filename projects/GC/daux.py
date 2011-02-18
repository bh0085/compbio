from Bio import SeqIO
import itertools as it
from matplotlib import pyplot as plt
from numpy import *

from compbio.biopython import biohilbert
from compbio.utils import memo as mem
from compbio.tests import drag_hilbert2 as dh2
import compbio.biopython.bhhash as bhhash
import compbio.config
import os

def gbk_elements(f):
  parsed = SeqIO.parse(f, 'genbank').next()



  nuclist = ['G','C','A','T']
  dinucs = []
  [[ dinucs.append( nuclist[j] + nuclist[i] ) for i in range(4)] for j in range(4)]
  dndict = {}.fromkeys(dinucs,0)
   
  
  fdict = {}
  fdicts = []
  ct = 0
  for f in parsed.features:    
    ct += 1
    idstr = 'ID'+str(ct)
    start, finish = f.location.start.position, \
        f.location.end.position
    seq = parsed.seq[start:finish]

    this_dict = {'g':seq.count('G'),
                   'c':seq.count('C'),
                   'len':len(seq),
                 'f':f
                   }
    for k in dndict.keys():
      this_dict[k] = seq.count(k)
    fdicts.append(this_dict)

 

  gcs = array([ float((x['g'] + x['c'])) / x['len'] for x in fdicts])
  skew = array([ float((x['g'] - x['c'])) /max([ x['g'] + x['c'],1.0]) for x in fdicts])

  elts = [{'end':fdicts[i]['f'].location.end.position,
           'start':fdicts[i]['f'].location.start.position,
           'gc':gcs[i],
           'skew':skew[i],
           'type':fdicts[i]['f'].type}
          for i in range(len(fdicts)) ]

  return elts

def gbk_speeded(elt_lists, speedfun):
  el_out = []
  for elts in elt_lists[0:5]:

    el_out.append(dict(children = [dict(
            gc = e['gc'],
            skew = e['skew'],
            type = e['type'],
            start = e['start'],
            length = e['end'] - e['start'] +1,
            id = 'blank_id',
            speed = speedfun(e)) for e in elts]))

  return el_out


def mycoplasma_elements():

  #Acceptable kwargs:
  #skeleton:
  # skel_lvls
  # skel_width
  # 
  #hilbert:
  # angle
  # lvls
  #
  #hash:
  # hash_spacer
  #
  #lambda
  root = compbio.config.root
  gdir = os.path.join(root, 'genomes/Bacteria')

  gbks = []
  for root, dirs, files in os.walk(gdir):
    for f in files:
      if 'gbk' in f:
        gbks.append(os.path.join(root, f))
  elt_lists = []
  for g in gbks:
    print g
    elt_lists.append(gbk_elements(g))
  
  return elt_lists

def run(reset = 0, **kwargs):

  if reset > 1:
    e_lists = mycoplasma_elements()
    mem.write('lists', e_lists,register = 'not')
  else:
    e_lists , sxs= mem.read('lists', register = 'not')
    if not sxs: raise Exception()


  kwargs['speedfun'] = kwargs.get('speedfun', lambda x: x['type'] == 'CDS' and 1.1 or 1)
  kwargs['widthfun'] =  kwargs.get('widthfun',lambda x: x['type'] == 'CDS' and 10 or 2)
  kwargs['colorfun'] =  kwargs.get('colorfun',lambda x: x['skew'] > 0 and 'red' or 'blue')

  if reset > 0: 
    e_lists2 = gbk_speeded(e_lists, kwargs['speedfun'])
    mem.write('lists2', e_lists2,register = 'el2', hardcopy = False)
  else:
    e_lists2 , sxs= mem.read('lists2',register = 'el2', hardcopy = False)
    if not sxs: raise Exception()        
    
  bh = biohilbert.BioHilbert(**kwargs)
  bh.init(e_lists2,  **kwargs)
  actor_kwargs = {
    'skeleton':kwargs.get('skel_lvls',2),
    'skelwidth':kwargs.get('skel_width',4),
    'qdraw':True
    }


  dh2.show_actor(bh, actor_kwargs = actor_kwargs) 

  

  return


