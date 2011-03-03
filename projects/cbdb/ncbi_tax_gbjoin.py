#!/usr/bin/env python

from numpy import *
import compbio.projects.cbdb.adapters as adapters
from sqlalchemy import Column, Integer, String, Unicode, ForeignKey
import compbio.config as config
import os

def get_maps():
  return dict(TaxGBJoin = {'gbid': 0,
                           'taxid':1})
                      
              
def get_tables():
  return [dict(name = 'TaxGBJoin',
               attrs = { 'gbid':Column(Integer,primary_key= True),
                         'taxid':Column(Integer)})]

def file_default():
  return config.dataPath('::ncbi/gi_taxid_nucl.dmp')


def splitfile():
  fopen = open(file_default())
  splits = 10000
  size = os.path.getsize(file_default())
  appx_ptrs = arange(0,splits,1,float) /splits * size
  real_ptrs = []
  for a in appx_ptrs:
    if a == 0:
      real_ptrs.append(int(a))
      continue
    fopen.seek(int(a))
    fopen.readline()
    real_ptrs.append(fopen.tell())
  starts = real_ptrs
  finishes = roll(list(real_ptrs),-1)
  finishes[-1] = size
  
  return starts, finishes


def fill_db(start = None, finish = None):
  filename = file_default()
  # a bit of a hack... there is only one table so it is just
  # the input filename
  table_defs = get_tables()
  column_maps =get_maps()
  table_names = map(lambda x: x['name'],table_defs)
  table_files = [filename]

  col_sep = '\t'
  rec_sep = '\n'

  adapters.fileMap2DB2('tax_gbs',
                      table_files,
                      table_names,
                      table_defs,
                      column_maps,
                      col_sep,
                      rec_sep,
                      321166*1,
                      reset = True)
                      
if __name__ == '__main__':
  print 'filling db'
  fill_db()
  exit(0)
