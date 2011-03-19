import compbio.projects.cbdb as cbdb
def phyla():
  dbi = cbdb.getName('taxdmp')
  phyla = dbi.S.q(dbi.Node).\
      filter_by(rank = 'phylum').all()
  return [sciname(p) for  p in phyla]

def sciname(node):
  sn= [name.name_txt 
          for name in node.names 
          if name.name_class == 'scientific name']
  return sn[0] if sn else ''
