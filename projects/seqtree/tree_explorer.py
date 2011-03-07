import compbio.projects.seqtree as sqt

def terminal_maps(tree, reset = False):
  gba_dbi, tax_dbi, tgb_dbi = sqt.gbaDBI(), sqt.taxDBI(), sqt.tgbDBI()
  tmaps = {}
  for t in tree.get_terminals():
    acc = sqt.clade_gbacc(t)
    gbid = gba_dbi.Session.query(gba_dbi.GBAcc).filter_by(accession = acc).all()
    if not gbid: 
      print '{0:20} for nodename: {1}'.format('failure', t.name)
      tmaps[t.name] = None
      continue
    print '{0:20} for nodename: {1}'.format('sxs', t.name)

    node = ncbiNodeForGBID(gbid[0])
    sciname = [x for x in node.names if x.name_class == 'scientific name'][0].name_txt
    tmaps[t.name] = node.id
  
  return tmaps
    
def ncbiNodeForGBID(gbid):
  gba_dbi, tax_dbi, tgb_dbi = sqt.gbaDBI(), sqt.taxDBI(), sqt.tgbDBI()
  taxid =  tgb_dbi.Session.query(tgb_dbi.TaxGBJoin).filter_by(gbid = gbid.id).one().taxid
  node = tax_dbi.Session.query(tax_dbi.Node).filter_by(id = taxid).one()    

  return node
