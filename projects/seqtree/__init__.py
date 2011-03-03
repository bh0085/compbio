from compbio.projects import SeqDB as sqd
from compbio import config
from compbio.utils import memo as mem, gbid_lookup as gbl
from Bio import Phylo, Align, AlignIO
import re
from numpy import * 
import compbio.projects.cbdb as cbdb
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


def rna4gbid( gbid, dbname = '16s'): 
  #print 'giving a random RNA because the taxonomy atabase is not yet created!'
  dbi = cbdb.getName(dbname)
  seq_num=  floor(1000 * random.random())
  seq = dbi.Session.query(dbi.Sequence)[seq_num].sequence
  return seq

def get_leaf_16s(clade):
  leaves = clade.get_terminals()
  rrnas = []
  random.seed(0)
  for l in leaves:
    gbacc= clade_gbacc(l)
    gbid = gbl.search_sorted(gbl.prefix(gbacc), gbacc)
    rrna = rna4gbid(gbid, dbname = '16s')
    rrnas.append(list(map(lambda x: ord(x),rrna)))

  arr = array([list(x) for x in rrnas])
  letters = sum( not_equal(arr, ord('-')), 0)

  ungapped_arr = arr[:,nonzero(letters)[0]]
  seq_letters = [''.join([chr(x) for x in y]).replace('-','-') for y in ungapped_arr]
  
  #there are about a thousand nonzero elements and really quite few 
  #gaps in the sequence that we get out of this method.

  align = Align.Generic.Alignment(Gapped(IUPAC.unambiguous_dna, "-"))
  for i in range(len(seq_letters)): align.add_sequence('SEQ%i  '%(i), seq_letters[i])
  AlignIO.write(align,open(config.dataPath('alignments/halo_16s.phylip'),'w'), 'phylip')
  AlignIO.write(align,open(config.dataPath('alignments/halo_16s.fasta'),'w'), 'fasta')
  AlignIO.write(align,open(config.dataPath('alignments/halo_16s.nexus'),'w'), 'nexus')
  Phylo.write(clade,open(config.dataPath('trees/halo_16s.newick'), 'w'), 'newick')

  raise Exception()
  print letters
  return rrnas
