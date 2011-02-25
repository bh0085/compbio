from sqlalchemy import Column, Integer, String, Unicode, ForeignKey, UniqueConstraint
from sqlalchemy.orm import relation
from compbio.projects import cbdb
from compbio import config
import os
import numpy as np
from Bio import SeqIO

class gbalign():
  def get_maps(self):
    return dict()

  def get_tables(self):
    #by convention, store alignment files in dataPath/alignments.
    return [dict(name = 'Alignment',
                 attrs={'id':Column(Integer, primary_key = True)}),
            dict(name = 'Sequence',
                 attrs={'id':Column(Integer, primary_key = True),
                        'name':Column(String),
                        'file_name':Column(String),
                        'file_offset':Column(Integer),
                        'sequence':Column(String),
                        'source_taxon':Column(Integer,index = True),
                        'source_organism':Column(String),
                        'gb_accession':Column(String),
                        'annotations':Column(String),
                        'alignmentid':Column(Integer, ForeignKey('alignment.id')),
                        'alignment':relation("Alignment",
                                             primaryjoin="Alignment.id==Sequence.alignmentid")
                        })
            ]
  def fill_db(self, name = '16s', reset = True):
    dbi = cbdb.getName(name, tables = self.get_tables(), 
                       reset = np.mod(reset,2))
    paths = []
    for r,ds, fs in os.walk(config.dataPath(os.path.join('alignments',name))):
      for f in fs:
        if 'gbk' in f: paths.append(os.path.join(r, f))
        
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
    
                           
        
