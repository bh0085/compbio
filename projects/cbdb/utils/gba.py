'''
  This contains a bunch of utilities useful for filling already
  created databases with additional metadata.

  The idea would be to create my databases with the default 
  scripts and then run the scripts in this file once my DB
  ecosystem is populated to throw in metadata that can only
  be set once all of the relevant databases have been created.
'''
import compbio.projects.cbdb as cbdb
def gba_fill_data(dbname):
  gba_map_gbids(dbname)
  gba_map_taxa(dbname)

def gba_map_gbids(dbname):
  dbi = cbdb.getName(dbname, tables =gb_alignment.get_tables())
  count = 0 
  fail_count =0
  sxs_count =0


  kind = 'dbsearch'
  if kind == 'filesearch':

    for s in dbi.Session.query(dbi.Sequence):    
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
        dbi.Session.flush()
        dbi.Session.commit()
        count = 0
      print prefix
  else:
    gbacc_dbi = cbdb.getName('gb_acc_idjoin' , gb_accid.get_tables())
    for s in dbi.Session.query(dbi.Sequence).all():
      try:
        gbid =  gbacc_dbi.Session.query(gbacc_dbi.GBAcc).\
            filter_by(accession = s.gb_accession).one().gbid
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
  dbi.Session.commit()

  print fail_count, sxs_count
def gba_map_taxa(dbname):
  seq_dbi = cbdb.getName(dbname,
                         tables = gb_alignment.get_tables())
  tgb_dbi = cbdb.getName('tax_gbs',
                         ncbi_tax_gbjoin.get_tables())

  count = 0 
  for  s in seq_dbi.Session.query(seq_dbi.Sequence):
    gbid = s.gb_id
    count += 1
    if gbid:
      taxid = tgb_dbi.Session.query(tgb_dbi.TaxGBJoin).\
          filter_by(gbid = gbid).one()
      s.source_taxon = taxid.taxid
    if count > 100:
      
      seq_dbi.Session.commit()
      print 'committing'
      count = 0
  seq_dbi.Session.commit()
