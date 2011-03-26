import compbio.projects.cbdb as cbdb
import compbio.utils.memo as mem

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
  if node == None: return None
  sn= [name.name_txt 
          for name in node.names 
          if name.name_class == 'scientific name']
  return sn[0] if sn else None

def taxon_with_name( rank, name, **kwargs ):
  def set_taxon_with_name(name = None, rank = None, **kwargs):
    assert name != None and rank != None
    all_p = get_rank(rank)
    p_node = [p for p in all_p if sciname(p) == name]
    assert len(p_node) == 1, 'Ambiguous phylum match?'
    p_node = p_node[0]
    return p_node
  return mem.getOrSet(set_taxon_with_name,
                      **dict(rank = rank,
                             name = name,
                             register = rank+'_'+name,
                             hardcopy = False,
                             on_fail = 'compute',
                             **kwargs))
