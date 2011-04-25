import compbio.projects.cbdb as cbdb

def get_rank(rankname = 'phylum'):
  dbi = cbdb.getName('taxdmp')
  rank = dbi.S.q(dbi.Node).\
      filter_by(rank = rankname).all()
  return rank

def get_root():
  dbi = cbdb.getName('taxdmp')
  return dbi.S.q(dbi.Node).\
      filter_by(id = 1).one()

def get_node(nodeid, dbi = None):
  if dbi == None: dbi = cbdb.getName('taxdmp')
  return dbi.S.q(dbi.Node).filter_by(id=nodeid).scalar()
  

def get_taxon(node, rank):
  if not node: return None
  rank_ancestor = [anc for anc in get_ancestors(node) if anc.rank == rank]
  return rank_ancestor[0] if rank_ancestor else None

def get_ancestors(node):
  l0 = [node]
  root = get_root()
  while l0[-1] != root:
    l0.append(l0[-1].parent)
    if len(l0) > 50: raise Exception()
  return l0

def sciname(node):
  sn= [name.name_txt 
          for name in node.names 
          if name.name_class == 'scientific name']
  return sn[0] if sn else ''
