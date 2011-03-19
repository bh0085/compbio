def map_paml(clade, dbi):
  '''
  Given a clade, map finds all of the accessions in the db
  that map to the terminals of the clade.

  For each mapped accesion, tracks down the appropriate
  sequence and writes newick and phylip files suitable for 
  ancestor inference with PAML.
  
  usage: map_paml(clade, dbname)
  '''

  seq_lists = []
  for t in clade.get_terminals():
    nodeid = t.m['taxnode'].id
    dbelts = dbi.Session.query(dbi.Sequence).\
        filter_by(source_taxon = nodeid).all()
    seq_lists.append([e.sequence for e in dbelts])
    
  taxids = [x [0] for x in dbi.S.q(dbi.Sequence.source_taxon).all()]
  tidxs = [t.m['taxnode'].id for t in clade.get_terminals()]
  matches = [t for t in taxids if t in tidxs]

  if len(matches) == 0:
    print '''No matches... did you remember to populate the
Database ID/Taxon fields using cbdb/utils/gba?'''
  raise Exception()

  



