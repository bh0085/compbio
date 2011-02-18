from Bio import AlignIO, SeqIO
from Bio.Align import MultipleSeqAlignment
from compbio import config
import os
from numpy import *
import pickle
import itertools as it

from compbio.utils import memo as mem

def index_gbks(path = os.path.join(config.root,'sequences/16s')):
  '''index a directory full of genbank files,
  pickling 'filename'.pindex for easy data retrieval.'''
  gbks = []
  basenames = []
  for r, d, files in os.walk(path):
    for f in files:
      if 'gbk' in f:
        gbks.append(os.path.join(r,f))
        basenames.append(os.path.splitext(f)[0])

  for i in range(len(gbks)):
    g = gbks[i]; b = basenames[i]
    subd = os.path.join(os.path.dirname(g),b)
    if not os.path.isdir(subd):
      os.mkdir(subd)
      
    f = open(g)
    prs = SeqIO.parse(f, 'genbank')
    multi = MultipleSeqAlignment([])
    starts, annotations, taxonomies, organisms ,gbids= [],[],[],[],[]
    print g,' Dir: ' ,subd
    while 1:

      this_start = f.tell() 
      try:
        seq = prs.next()
      except:
        break
      starts.append(this_start)
      #annotations.append(seq.annotations)

      #there appear to be some issues with taxonomy names.
      tax = map(lambda x: x.replace('"',''), seq.annotations['taxonomy'])
      taxonomies.append(tax)
      gbids.append(seq.annotations['comment'][9:])
      organisms.append(seq.annotations['organism'])


    n = len(starts)
    #pickle.dump([[annotations[i], starts[i]] for i in range(n)]
    #            , open(os.path.join(subd,'annotations'),'w'))

    #sort by gbids...
    tax_srt = argsort(map(lambda x: x.__repr__(), taxonomies))
    pickle.dump([[taxonomies[tax_srt[i]],starts[tax_srt[i]] ]for i in range(n)]
                , open(os.path.join(subd,'taxonomies'),'w'))

    #sort by gbids...
    gbsrt = argsort(gbids)
    gbsrted = array(gbids)[gbsrt]
    gbstarts = array(starts)[gbsrt]
    pickle.dump([array(gbsrted), array(gbstarts)], 
                open(os.path.join(subd,'gbids'),'w'))

def make_tax_idxs(path = os.path.join(config.root,'sequences/16s'),reset = False):
    all_idxs = {}
    for root, dirs, files in os.walk(path):
      if 'taxonomies' in files and 'gbids' in files:
          t = pickle.load(open(os.path.join(root, 'taxonomies')))
          g = pickle.load(open(os.path.join(root, 'gbids')))
                 
          gstart_sort = argsort(g[1])
          gstarts = g[1][gstart_sort]
          
          phyla = {}
          #sort by phyla
          
          taxons = {1:'domains',
                  2:'phyla',
                  3:'classes',
                  4:'orders',
                  5:'families',
                  6:'genii',
                  7:'species'}
          for i, taxon in taxons.iteritems():
            print root, taxon
            for k, grp in it.groupby(t, lambda x: len(x[0])<= i and 'NA' or x[0][i]):
              gbs = []
              for elt in grp:
                start = elt[1]
                gidx = gstart_sort[searchsorted(gstarts,start)]
                gbs.append( g[0][gidx] )
              phyla[k] = gbs
            
            pickle.dump(gbs, open(os.path.join(root, taxon),'w'))
          

def get_gb_idxs(path = os.path.join(config.root,'sequences/16s'),reset = False):
  if not reset:
    out, sxs = mem.read('idxs')
    assert sxs
  else:
    all_idxs = {}
    for root, dirs, files in os.walk(path):
      for f in files:
        if f == 'gbids':
          all_idxs[os.path.join(root, f)] = \
              pickle.load(open(os.path.join(root, f)))
    mem.write('idxs', all_idxs)
    out = all_idxs
  return out

  
def get_tax_idxs(path = os.path.join(config.root, 'sequences/16s'), reset = False):
  if not reset:
    out, sxs = mem.read('idxs')
    assert sxs
  else:    
    all_idxs = {}
    for root, dirs, files in os.walk(path):
      for f in files:
        if f == 'taxonomies':
          all_idxs[os.path.join(root, f)] = \
              pickle.load(open(os.path.join(root, f)))
    mem.write('idxs', all_idxs)
    out = all_idxs
  return out

def gbid_record(gbid, reset = False):
  idxs = get_gb_idxs( reset = mod(reset, 2))
  
  file_found = None
  file_start = -1
  for k, vals in idxs.iteritems():
    match = searchsorted(vals[0], gbid)
    if vals[0][match] == gbid:
      file_found = k
      file_start = vals[1][match]
      break
  assert file_found
  f = open(os.path.dirname(file_found) + '.gbk')
  f.seek(file_start)
  rec = SeqIO.parse(f, 'genbank').next()
  return rec           
    
def tax_records(tax_str,max_hits = 100, reset = False):
  tax_idxs = get_tax_idxs(reset = mod(reset,2))
  for k , v  in tax_idxs:
    for elt in v:
      if tax_str in elt.__repr__():
        print elt
      

    break
  
  
