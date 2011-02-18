from Bio import AlignIO
from compbio import config
import os

def index_gbks(path = os.path.join(config.root,'sequences/16s')):
  '''index a directory full of genbank files,
  pickling 'filename'.pindex for easy data retrieval.'''
  gbks = []
  for r, d, files in os.walk(path):
    for f in files:
      if 'gbk' in f:
        gbks.append(f)

  
  print gbks
