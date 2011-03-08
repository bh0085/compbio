from compbio.projects.cbdb import *
import compbio.projects.cbdb as cbdb
def fill_data(dbname):
  map_gbids(dbname)
  map_taxa(dbname)

def map_gbids(dbname):
  dbi = cbdb.getName(dbname,
                     tables =gb_alignment.get_tables())
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
        dbi.Session.commit()
        count = 0
      print prefix
  else:
    gbacc_dbi = cbdb.getName('gb_accjoin' , gb_accid.get_tables())
    raise Exception()

  print fail_count, sxs_count
def map_taxa(dbname):
  seq_dbi = cbdb.getName(dbname,
                         tables = gb_alignment.get_tables())
  tax_dbi = cbdb.getName('taxdmp',
                         ncbi_tax.get_tables())

  for  s in seq_dbi.Session.query(seqdbi.Sequence):
    gbid = s.gb_id
    raise Exception()
