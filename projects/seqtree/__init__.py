from compbio.projects import SeqDB as sqd
from compbio import config
from compbio.utils import memo as mem, gbid_lookup as gbl
from Bio import Phylo, Align, AlignIO
import re, itertools as it
from numpy import * 
import compbio.projects.cbdb as cbdb
from Bio.Alphabet import IUPAC, Gapped
import numpy as np


def init(**kwargs):
  '''
  Read in the 16s tree of life and a random clade corresponding to the
  halobacteria.

  At each node, sets metadata from the databases that I have grabbed.
  Metadata (node.m) for terminal nodes includes:
    taxnode   -- ncbi taxon of the node
    gbacc     -- genbank accession number of the 16s for the node
    gbid      -- genbank id of the 16s for the node
    
  inputs:
    reset [False]

  output:
    tree  <biopython tree>, the entire 16s tree of life
    halo  <biopython clade>, a clade of the tree of life

  usage:
    tree, halo = init()
'''
  
  print 'testing...'
  def setTree(**kwargs):
    nwk = Phylo.read(config.dataPath('sequences/16s.newick'),"newick")
    for n in it.chain(nwk.get_terminals(),nwk.get_nonterminals()): n.m = {}
    db_metadata(nwk)
    print "SETTING TREE!!!"
    return nwk
  
  return mem.getOrSet(setTree,
                      **mem.rc( kwargs, 
                                name = kwargs.get('name', 'default_tree'),
                                on_fail = 'compute',
                                register = 'init'))
  
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


def db_metadata(clade):
  '''
  Get an ncbi genealogy for a clade - e.g: the minimal ncbi node containing
  every terminal of a clade as well as the ncbi nodes at each leaf.

  inputs:
    clade: <biopython clade>
    
  outputs:
    genealogy: the shared ncbi genealogy of every terminal
    nodes:     ncbi taxnodes for every terminal
'''
  tax_dbi = cbdb.getName('taxdmp')
  tgb_dbi = cbdb.getName('tax_gbs')
  gba_dbi = cbdb.getName('gb_acc_idjoin')
  
  #terminal_gbaccs = map(lambda x: clade_gbacc(x), clade.get_terminals())  

  for idx, t in enumerate(clade.get_terminals()):
    try:
      t.m['id'] = idx
      t.m['gbacc'] = clade_gbacc(t)
      t.m['gbid'] = gba_dbi.Session.query(gba_dbi.GBAcc).\
          filter_by(accession = t.m['gbacc']).one().gbid
      taxid = tgb_dbi.Session.query(tgb_dbi.TaxGBJoin).\
          filter_by(gbid = t.m['gbid']).one().taxid
      t.m['taxid'] = taxid
    except: pass
  max_idx = idx
  for idx, t in enumerate(clade.get_nonterminals()):
    t.m['id'] = max_idx + idx


def taxRoot():
  tax_dbi = cbdb.getName('taxdmp')
  root_node = tax_dbi.Session.query(tax_dbi.Name).filter_by(name_txt  = 'root').one().node
  return root_node

def getGenealogy(node):
  root_node =taxRoot()
  dbi = cbdb.getName('taxdmp')

  path = []
  cur = node
  while cur != root_node:
    path.append(cur)
    cur = cur.parent
  path.append(cur)
  return path[::-1]

def rna4gbid( gbid, dbname = '16s'): 
  #print 'giving a random RNA because the taxonomy atabase is not yet created!'
  dbi = cbdb.getName(dbname)
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

  raise Exception()
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

  #t0 = Phylo.BaseTree.Tree(p0)
  Phylo.write(t0,open(config.dataPath('trees/halo_16s.newick'), 'w'), 'newick')
