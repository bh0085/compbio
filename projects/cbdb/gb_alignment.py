from sqlalchemy import Column, Integer, String, Unicode, ForeignKey, UniqueConstraint
from sqlalchemy.orm import relation
from compbio.projects import cbdb
from compbio import config
import os
import numpy as np
from Bio import SeqIO

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
                      'gb_accession':Column(String),
                      'gb_accession_range':Column(String),
                      'gb_id':Column(String),
                      'annotations':Column(String),
                      'alignmentid':Column(Integer, ForeignKey('alignment.id')),
                      'alignment':relation("Alignment",
                                           primaryjoin="Alignment.id==Sequence.alignmentid")
                      })
          ]

def fill_from_rfam_stk( p, reset = True):

  aname = os.path.basename(p)
  dbi = cbdb.getName(aname, tables = get_tables(), 
                     reset = np.mod(reset,2))
      
  fopen = open(p)
  a = dbi.Alignment(file_name = aname)
  dbi.Session.add(a)
  dbi.Session.commit()

  count = 0
  for rec in SeqIO.parse(fopen, 'stockholm'):
    source_taxon = None
    source_organism = None

    acc = rec.annotations['accession']
    accid, accrange = acc.split('/')

    seq = dbi.Sequence(name = rec.name,
                       file_name = p,
                       file_offset = fopen.tell(),
                       sequence = rec.seq.__str__(),
                       source_taxon = source_taxon,
                       source_organism = source_organism,
                       gb_accession = accid,
                       gb_accession_range = accrange,
                       gb_id = None,
                       annotations = rec.annotations.__str__(),
                       alignmentid = 0)

    dbi.Session.add(seq)
    if np.mod(count, 100) == 0:
      print count, p , seq.source_organism
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
    
                           
        
