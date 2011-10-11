import cb.config as cfg
import numpy as np
from numpy import *
import os
import cb.utils.memo as mem

import Bio.SeqIO as sio
import cb.external.GFF as GFF

def peaks(**kwargs):
    def setPeaks(**kwargs):
        root = cfg.dataPath('wormchip')
        files = [os.path.join(root, f) for f in os.listdir(root)]
        parser = GFF.GFFParser()
        datum =list( parser.parse(files))
        
        out = dict([(os.path.basename(files[i]),
                     [{'start':f.location.start.position,
                           'end':f.location.end.position,
                           'type':f.type,
                           'score':f.qualifiers['score'][0],
                           'qValue':f.qualifiers['qValue'][0]}
                      for f in d.features])
                    for i,d in enumerate(datum)])
        
        return out
    return mem.getOrSet(setPeaks, 
                        **mem.rc(kwargs,
                                 hardcopy = True)
                        )

def chromosome_names():
    return  ['CHR_I', 'CHR_II', 'CHR_III'
             ,'CHR_IV','CHR_V','CHR_X']
def chromosome_offsets(**kwargs):
    def set_chromosome_offsets(**kwargs):

      lens =[]
      names = chromosome_names()
      for name in names:
         root = cfg.dataPath('/data/genomes/Caenorhabditis_elegans')
         fdir = os.path.join(root,name)
         for r, d, files in os.walk(fdir):
             for f in files:
                 if '.gb' in f:
                     fopen = open(os.path.join(r,f))
                     break
      
         gb = list(sio.parse(fopen, 'genbank'))[0]
         fopen.close()
         lens.append( gb.features[0].location.end.position)

      offsets = {}
      cur_ofs = 0
      for i, l in enumerate(lens):
        offsets[names[i]] = cur_ofs
        cur_ofs += l 
      return offsets
    return mem.getOrSet(set_chromosome_offsets, 
                        **mem.rc(kwargs, hardcopy = True))
    
def parse_genes(**kwargs):

    def set_genes(**kwargs):

      lens =[]
      all_genes = {}
      names = chromosome_names()
      for name in names:
         root = cfg.dataPath('/data/genomes/Caenorhabditis_elegans')
         fdir = os.path.join(root,name)
         for r, d, files in os.walk(fdir):
             for f in files:
                 if '.gb' in f:
                     fopen = open(os.path.join(r,f))
                     break
      

         gb = list(sio.parse(fopen, 'genbank'))[0]
         genes = [f for f in gb.features if f.type== 'gene']
         all_genes[name] = genes
         fopen.close()
      return all_genes
    return mem.getOrSet(set_genes, 
                        **mem.rc(kwargs,hardcopy = True))

