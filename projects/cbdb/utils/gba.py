#!/usr/bin/env python
'''
  This contains a bunch of utilities useful for filling already
  created databases with additional metadata.

  The idea would be to create my databases with the default 
  scripts and then run the scripts in this file once my DB
  ecosystem is populated to throw in metadata that can only
  be set once all of the relevant databases have been created.
'''
from sqlalchemy import func
import compbio.projects.cbdb
def gba_fill_data(dbname):
  gba_map_gbids(dbname)
  gba_map_taxa(dbname)

def gba_map_gbids(dbname):
  cbdb = compbio.projects.cbdb

  dbi = cbdb.getName(dbname)
  count = 0 
  fail_count =0
  sxs_count =0


  gbacc_dbi = cbdb.getName('gb_acc_idjoin')
  slc_n = 10000
  slc_ofs = 0
  max_ofs =dbi.S.q(func.max(dbi.Sequence.id)).scalar() 
  while slc_ofs < max_ofs:
    for s in dbi.Session.query(dbi.Sequence).\
          filter(dbi.Sequence.id>=slc_ofs).\
          filter(dbi.Sequence.id<(slc_ofs+slc_n)):
      try:
        gbid =  gbacc_dbi.Session.query(gbacc_dbi.GBAcc).\
            filter_by(accession = s.gb_accession).first().gbid
        s.gb_id = gbid
        sxs_count += 1
      except:
        print 'failed!'
        fail_count +=1
      count += 1
      if count > 100:
        dbi.Session.commit()
        count = 0
        print 'adding'
    slc_ofs += slc_n

  dbi.Session.commit()

  print fail_count, sxs_count
def gba_map_taxa(dbname):
  cbdb = compbio.projects.cbdb
  dbi = cbdb.getName(dbname)
  tgb_dbi = cbdb.getName('tax_gbs')

  count = 0 
  slc_n = 10000
  slc_ofs = 0
  max_ofs =dbi.S.q(func.max(dbi.Sequence.id)).scalar() 
  while slc_ofs < max_ofs:
    for  s in dbi.Session.query(dbi.Sequence).\
          filter(dbi.Sequence.id>=slc_ofs).\
          filter(dbi.Sequence.id<(slc_ofs+slc_n)):
      gbid = s.gb_id
      count += 1
      if gbid:
        taxid = tgb_dbi.Session.query(tgb_dbi.TaxGBJoin).\
            filter_by(gbid = gbid).first()
        s.source_taxon = taxid.taxid
      if count > 100:
        
        dbi.Session.commit()
        print 'committing'
        count = 0
    slc_ofs += slc_n
  dbi.Session.commit()

import sys
if __name__ == '__main__':
  args = sys.argv[1:]
  dbname = args[0]
  print 'making metadata for {0}'.format(dbname)
  gba_fill_data(dbname)
  
