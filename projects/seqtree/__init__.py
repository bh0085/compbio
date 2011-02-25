from compbio.projects import SeqDB as sqd
from compbio import config
from compbio.utils import memo as mem
from Bio import Phylo

def init(reset = False):
  dbi = sqd.dbname('16s')
  
  if reset: 
    nwk = Phylo.read(config.absPath('sequences/16s.newick'),"newick")
    mem.write(value = nwk, register = '16s_tree')
  else: 
    nwk, sxs = mem.read(register = '16s_tree')
    assert sxs

  halo = Phylo.BaseTree.Tree(nwk.find_clades(dict(name = 'Halobacteria')).next())
  return halo
  
def clade_gbid(clade):
  gbid = re.compile('__([^_]*)__').search(clade.name)
  return gbid.group(1)
