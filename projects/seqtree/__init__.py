from compbio.projects import SeqDB as sqd
from compbio import config
from compbio.utils import memo as mem, gbid_lookup as gbl
from Bio import Phylo, Align, AlignIO
import re
from numpy import * 
import compbio.projects.cbdb as cbdb
from compbio.projects.cbdb import *
import compbio.projects.cbdb.gb_alignment as gba
from Bio.Alphabet import IUPAC, Gapped

def init(reset = False):
  
  if reset: 
    nwk = Phylo.read(config.dataPath('sequences/16s.newick'),"newick")
    mem.write(value = nwk, register = '16s_tree')
  else: 
    nwk, sxs = mem.read(register = '16s_tree')
    assert sxs

  halo = nwk.find_clades(dict(name = 'Halobacteria')).next()
  return nwk, halo
  
def clade_gbacc(clade):
  gbid = re.compile('__([^_]*)__').search(clade.name)
  return gbid.group(1)


def write_nwk(tree):
  #tree = Phylo.BaseTree.Tree(node)
  fout =config.dataPath('::trees/halo.newick')
  Phylo.write(tree,open(fout,'w'), 'newick')

def clade_taxa(clade):
  terminal_gbaccs = map(lambda x: clade_gbacc(x), clade.get_terminals())
  
  tax_dbi = cbdb.getName('taxdmp', ncbi_tax.get_tables())
  tgb_dbi = cbdb.getName('tax_gbs', ncbi_tax_gbjoin.get_tables())
  gba_dbi = cbdb.getName('gb_accjoin', gb_accid.get_tables())
  
  ids = []
  for t in terminal_gbaccs:
    gbid = gba_dbi.Session.query(gba_dbi.GBAcc).filter_by(accession = t).all()
    if not gbid:
      continue
    
    taxid = tgb_dbi.Session.query(tgb_dbi.TaxGBJoin).filter_by(gbid = gbid[0].id).one().taxid
    node = tax_dbi.Session.query(tax_dbi.Node).filter_by(id = taxid).one()    
    g = getGenealogy(node)
    print map(lambda x: x.names[map(lambda y: y.name_class, x.names).index('scientific name')].name_txt, g)    

    ids.append(gbid)

  raise Exception()

def taxRoot():
  tax_dbi = taxDBI()
  root_node = tax_dbi.Session.query(tax_dbi.Name).filter_by(name_txt  = 'root').one().node
  return root_node
  
def gbaDBI():
  gba_dbi = cbdb.getName('gb_accjoin', gb_accid.get_tables())
  return gba_dbi

def tgbDBI():
  tgb_dbi = cbdb.getName('tax_gbs', ncbi_tax_gbjoin.get_tables())
  return tgb_dbi

def taxDBI():
  tax_dbi = cbdb.getName('taxdmp', ncbi_tax.get_tables())
  return tax_dbi

def getGenealogy(node):
  root_node =taxRoot()
  dbi = taxDBI()

  path = []
  cur = node
  while cur != root_node:
    path.append(cur)
    cur = cur.parent
  path.append(cur)
  return path[::-1]

def rna4gbid( gbid, dbname = '16s'): 
  #print 'giving a random RNA because the taxonomy atabase is not yet created!'
  dbi = cbdb.getName(dbname, gba.get_tables())
  seq_num=  floor(1000 * random.random())
  seq = dbi.Session.query(dbi.Sequence)[seq_num].sequence
  return seq

def get_leaf_16s(clade):
  cltree = Phylo.BaseTree.Tree(clade)
  leaves = clade.get_terminals()
  
  
  l0 = leaves[0]
  p0 = clade.get_path(l0)[-3]
  siblings = p0.get_terminals()

  rrnas = []
  random.seed(5)

  t0 = Phylo.BaseTree.Tree(p0)
  
  ct = 0
  names = []
  for l in t0.get_terminals():
    gbacc= clade_gbacc(l)
    gbid = gbl.search_sorted(gbl.prefix(gbacc), gbacc)
    rrna = rna4gbid(gbid, dbname = '16s')
    rrnas.append(list(map(lambda x: ord(x),rrna)))
    l.name = 'SEQ%i  '%(ct)
    names.append(l.name)
    ct += 1
  arr = array([list(x) for x in rrnas])
  letters = sum( not_equal(arr, ord('-')), 0)

  ungapped_arr = arr[:,nonzero(letters)[0]]
  seq_letters = [''.join([chr(x) for x in y]).replace('-','-') for y in ungapped_arr]
  
  #there are about a thousand nonzero elements and really quite few 
  #gaps in the sequence that we get out of this method.

  align = Align.Generic.Alignment(Gapped(IUPAC.unambiguous_dna, "-"))
  for i in range(len(seq_letters)): align.add_sequence(names[i], seq_letters[i])
  AlignIO.write(align,open(config.dataPath('alignments/halo_16s.phylip'),'w'), 'phylip')
  AlignIO.write(align,open(config.dataPath('alignments/halo_16s.fasta'),'w'), 'fasta')
  AlignIO.write(align,open(config.dataPath('alignments/halo_16s.nexus'),'w'), 'nexus')

  t0 = Phylo.BaseTree.Tree(p0)
  Phylo.write(t0,open(config.dataPath('trees/halo_16s.newick'), 'w'), 'newick')

  print map(lambda x: x.seq.__str__(), align)
  raise Exception()
  print letters
  return rrnas
