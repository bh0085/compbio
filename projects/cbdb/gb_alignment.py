from sqlalchemy import Column, Integer, String, Unicode, ForeignKey, UniqueConstraint
from sqlalchemy.orm import relation
from compbio.projects import cbdb
import compbio.utils.gbid_lookup as gbl
from compbio import config
import os, re
import numpy as np
from Bio import SeqIO
import compbio.projects.cbdb.ncbi_tax as nct

def fill_data(dbname):
  map_gbids(dbname)
  map_tax(dbname)

def map_gbids(dbname):
  dbi = cbdb.getName(dbname,
                     tables =get_tables())
  count = 0 
  fail_count =0
  sxs_count =0
  for s in dbi.Session.query(dbi.Sequence):
    print s.gb_accession
    acc = s.gb_accession.split('.')[0]
    prefix = re.compile('[A-Z]*').search(acc).group()
    try:
      gbid =  gbl.search_sorted(prefix, acc)
      s.gb_id = gbid
      sxs_count +=1
    except Exception, e:
      fail_count +=1
    count += 1
    if count > 100:
      dbi.Session.commit()
      count = 0
      print prefix

def map_taxa(dbname):
  seq_dbi = cbdb.getName(dbname,
                     tables = get_tables())
  tax_dbi = cbdb.getName('taxdmp',
                         nct.get_tables())

  for  s in seq_dbi.Session.query(seqdbi.Sequence):
    gbid = s.gb_id
    raise Exception()

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
                      'gb_id':Column(Integer, index = True),
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
    acc = rec.annotations['accession']
    accid, accrange = acc.split('/')

    seq = dbi.Sequence(name = rec.name,
                       file_name = p,
                       file_offset = fopen.tell(),
                       sequence = rec.seq.__str__(),
                       gb_accession = accid,
                       gb_accession_range = accrange,
                       gb_id = None,
                       annotations = rec.annotations.__str__(),
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
      seq = dbi.Sequence(name = rec.name,
                         file_name = p,
                         file_offset = last_ofs,
                         sequence = rec.seq.__str__(),
                         gb_accession = rec.id,
                         annotations = rec.annotations.__str__(),
                         alignment = a,
                         source_taxon = src_taxon
                         )
      dbi.Session.add(seq)
      last_ofs = fopen.tell()
      if np.mod(count, 1000) == 0:
        print count, p, seq.source_organism
        dbi.Session.commit()
      count += 1
    dbi.Session.commit

        

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
    
                           
        
