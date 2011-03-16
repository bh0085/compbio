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
  '''
  Read in the 16s tree of life and a random clade corresponding to the
  halobacteria.

  inputs:
    reset [False]

  output:
    tree  <biopython tree>, the entire 16s tree of life
    halo  <biopython clade>, a clade of the tree of life

  usage:
    tree, halo = init()
'''
  
  if reset: 
    nwk = Phylo.read(config.dataPath('sequences/16s.newick'),"newick")
    mem.write(value = nwk, register = '16s_tree')
  else: 
    nwk, sxs = mem.read(register = '16s_tree')
    assert sxs

  halo = nwk.find_clades(dict(name = 'Halobacteria')).next()
  return nwk, halo
  
def clade_gbacc(clade):
  '''
  Run on a clade from the tree of life I have been using, returns a
  genbank accession number for any clade corresponding to a 16s accession.

  inputs:
    clade: <biopython clade>
   
  outputs:
    gbacc: genbank accession
'''
  gbid = re.compile('__([^_]*)__').search(clade.name)
  return gbid.group(1)


def clade_taxa(clade):
  '''
  Get an ncbi genealogy for a clade - e.g: the minimal ncbi node containing
  every terminal of a clade as well as the ncbi nodes at each leaf.

  inputs:
    clade: <biopython clade>
    
  outputs:
    genealogy: the shared ncbi genealogy of every terminal
    nodes:     ncbi taxnodes for every terminal
'''
  terminal_gbaccs = map(lambda x: clade_gbacc(x), clade.get_terminals())
  
  tax_dbi = cbdb.getName('taxdmp', ncbi_tax.get_tables())
  tgb_dbi = cbdb.getName('tax_gbs', ncbi_tax_gbjoin.get_tables())
  gba_dbi = cbdb.getName('gb_accjoin', gb_acc_idjoin.get_tables())
  
  ids = []
  taxnodes = []
  genealogies = []
  for t in terminal_gbaccs:
    gbid = gba_dbi.Session.query(gba_dbi.GBAcc).filter_by(accession = t).one()
    if not gbid:
      raise Exception('uh oh... is the gbacc db still borked?')
    
    taxid = tgb_dbi.Session.query(tgb_dbi.TaxGBJoin).filter_by(gbid = gbid.id).one().taxid
    node = tax_dbi.Session.query(tax_dbi.Node).filter_by(id = taxid).one()    
    g = getGenealogy(node)
    print map(lambda x: x.names[map(lambda y: y.name_class, x.names).index('scientific name')].name_txt, g)    

    ids.append(gbid)
    taxnodes.append(node)
    genealogies.append(g)

  mlen = np.min([len(g) for g in genealogies])
  shared = 0
  for i in range(mlen):
    if len(nonzero([g[i] != genealogies[0][i] for g in genealogies])[0]):
      break
  fully_shared = i

  return(genealogies[0:fully_shared], taxnodes)

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
