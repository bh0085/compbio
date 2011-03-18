#!/usr/bin/env python

from sqlalchemy import Column, Integer, String, Unicode, ForeignKey, UniqueConstraint
from sqlalchemy.orm import relation
import compbio.projects.cbdb
import compbio.utils.gbid_lookup as gbl
from compbio import config
import os, re, sys
import numpy as np
from Bio import SeqIO
import compbio.utils.pbar as pbar
import compbio.projects.cbdb
import simplejson as sjson


def get_maps():
  return dict()

def get_tables():
  #by convention, store alignment files in dataPath/alignments.
  return [dict(name = 'Alignment',
               attrs={'id':Column(Integer, primary_key = True),
                      'file_name':Column(String)}),
          dict(name = 'Sequence',
               attrs={'id':Column(Integer, primary_key = True),
                      'file_name':Column(String),
                      'file_offset':Column(String),
                      'name':Column(String),
                      'sequence':Column(String),
                      'source_taxon':Column(Integer,index = True),
                      'source_organism':Column(String),
                      'gb_accession':Column(String,index = True),
                      'gb_accession_version':Column(Integer),
                      'gb_accession_range':Column(String),
                      'gb_id':Column(Integer, index = True),
                      'annotations':Column(String),
                      'alignmentid':Column(Integer, ForeignKey('alignment.id')),
                      'alignment':relation("Alignment",
                                           primaryjoin="Alignment.id==Sequence.alignmentid")
                      })
          ]

def fill_from_rfam_stk( p, reset = True):
  cbdb = compbio.projects.cbdb
  aname = os.path.basename(p)
  dbi = cbdb.getName(aname, tables = get_tables(), 
                     reset = np.mod(reset,2))
      
  fopen = open(p)
  a = dbi.Alignment(file_name = aname)
  dbi.Session.add(a)
  dbi.Session.commit()

  count = 0
  for rec in SeqIO.parse(fopen, 'stockholm'):
    acc = rec.annotations['accession']
    accidv, accrange = acc.split('/')
    acv_split =  accidv.split('.')
    accid = acv_split[0]
    accid_version = (lambda x: len(x) == 1 and 1 or x[1])(acv_split)
    
    ann = sjson.dumps(rec.annotations, default = lambda x: x.__str__())
    seq = dbi.Sequence(name = rec.name,
                       file_name = p,
                       file_offset = fopen.tell(),
                       sequence = rec.seq.__str__(),
                       gb_accession = accid,
                       gb_accession_version = accid_version,
                       gb_accession_range = accrange,
                       gb_id = None,
                       annotations = ann,
                       alignment = a
                       )
    

    dbi.Session.add(seq)
    if np.mod(count, 100) == 0:
      print count, p , seq.source_organism
      dbi.Session.commit()
    count += 1
  dbi.Session.commit()
  
                         
def fill_all_rdb16s(reset = True):
  paths = []
  for r, ds, fs in os.walk(config.dataPath('alignments/16s')):
    for f in fs:
      if '.gbk' in f:
        paths.append(os.path.join(r,f))
  cbdb = compbio.projects.cbdb
  dbi = cbdb.getName('16s',
                     tables = get_tables(),
                     reset = np.mod(reset, 2))
  last_ofs = 0
  for p in paths:
    fopen = open(p)
    a = dbi.Alignment(file_name =config.dataURL(p))
    dbi.Session.add(a)
    dbi.Session.commit()
    count = 0 

    for rec in SeqIO.parse(fopen, 'genbank'):
      try:
        src_taxon = rec.features[0].qualifiers['db_xref'][0][6:]
      except Exception, e:
        src_taxon = None

      ann = sjson.dumps(rec.annotations, default = lambda x: x.__str__())
      seq = dbi.Sequence(name = rec.name,
                         file_name = p,
                         file_offset = last_ofs,
                         sequence = rec.seq.__str__(),
                         gb_accession = rec.id,
                         gb_accession_version = 1,
                         gb_id = None,
                         annotations = ann,
                         alignment = a,
                         source_taxon = src_taxon
                         )
      dbi.Session.add(seq)
      last_ofs = fopen.tell()
      if np.mod(count, 1000) == 0:
        print count, p, seq.source_organism
        dbi.Session.commit()
      count += 1
    dbi.Session.commit()

        

def fill_db( name = '16s', reset = True):
  dbi = cbdb.getName(name, tables = get_tables(), 
                     reset = np.mod(reset,2))
      
  count = 0 
  for p in paths:
    fopen = open(p)
    a = dbi.Alignment()
    dbi.Session.add(a)
    dbi.Session.commit()
    for rec in SeqIO.parse(fopen, 'genbank'):
      f0 = rec.features[0]
      if f0.type == 'source':
        source_taxon = f0.qualifiers['db_xref'][0][6:]
        source_organism=f0.qualifiers['organism'][0]
      else:
        source_taxon = None
        source_organism = None

      seq = dbi.Sequence(name = rec.name, 
                         file_name = p,
                         file_offset = fopen.tell(),
                         sequence = rec.seq.__str__(),
                         source_taxon = source_taxon,
                         source_organism = source_organism,
                         gb_accession = rec.id,
                         annotations = rec.annotations.__str__(),
                         alignmentid = 0)

      dbi.Session.add(seq)
      if np.mod(count, 1000) == 0:
        print count, p , seq.source_organism
        dbi.Session.commit()
      count += 1
    dbi.Session.commit()
    
                           
        
if __name__ == '__main__':
  args = sys.argv[1:]
  if args[0] == '16s':
    fill_all_rdb16s(reset = True)
