import compbio.projects.cbdb
import compbio.utils.memo as mem
import compbio.projects.cbdb.helpers.ncbi as ncbi

def get_seqs(dbname, **kwargs):
  def set_seqs(**kwargs):
    cbdb = compbio.projects.cbdb
    dbname = kwargs['dbname']
    dbi = cbdb.getName(dbname)
    nodes = dbi.S.q(dbi.Sequence).all()
    return nodes
  kwnew =  mem.rc(kwargs,hardcopy = False,
                  name = dbname, on_fail = 'compute',
                  dbname = dbname)
  return mem.getOrSet(set_seqs, **kwnew)

def get_taxon_forsome(nodes,rank,set_name = 'default_setname',**kwargs):
  def set_taxon_forsome(nodes = None, rank = None,**kwargs):
    assert nodes != None and rank != None
    taxon = [ncbi.get_taxon(node, rank = rank)
             if node else None for node in nodes]
    return taxon
  
  return mem.getOrSet(set_taxon_forsome,
                      **mem.rc(kwargs,
                               nodes = nodes,
                               rank = rank,
                               on_fail = 'compute',
                               hardcopy = False,
                               register= set_name + rank))
  

def get_taxon_forall(aliname,
                     rank = None, 
                     **kwargs):
  def setTaxon(aliname = None, rank = None,**kwargs):
    assert aliname != None and rank != None
    nodes = get_taxnodes(aliname,**mem.skw(**kwargs))
    taxon = [ncbi.get_taxon(node, rank=rank) 
             if node else None for node in nodes]
    return taxon
  return mem.getOrSet(setTaxon,
                      **mem.rc(kwargs,
                              aliname = aliname,
                              rank = rank,
                              on_fail = 'compute',
                              hardcopy = False,
                              register = aliname + rank))

def get_taxnodes(dbname, **kwargs):
  def set_taxnodes(**kwargs):
    
    all_seqs = get_seqs(dbname,**mem.skw(**kwargs))
    seq_taxa = [s.source_taxon 
                   if s.source_taxon else None 
                   for s in all_seqs]
    alinodes = [ncbi.get_node(s) if s != None else None for s in seq_taxa]
    return alinodes
  return mem.getOrSet(set_taxnodes,
                      **mem.rc(kwargs,
                               on_fail = 'compute',
                               hardcopy = False, 
                               register = dbname))
