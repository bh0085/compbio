from Bio import Phylo, SeqIO, AlignIO
from compbio import config
import os
import re

def init():


  nwk = get_16s_nwk()
  ali = get_16s_ali()
  genomes = get_myc_genomes()
  return nwk, ali, genomes


def get_16s_nwk():
  nwk = Phylo.read(os.path.join(config.root,
                               'sequences/16s.newick') ,"newick")
  return nwk

def get_16s_ali():
  f = os.path.join(config.root, 'sequences/rdb_13s_archaea.gen')
  seqs = SeqIO.parse(f,'genbank')
  als = AlignIO.parse(f,'genbank').next()

  return als

def clade_ancestors():
  pass

def clade_gbid(clade):
  gbid = re.compile('__([^_]*)__').search(clade.name)
  return gbid.group(1)

def get_16s_gbid(gbk):
  return gbk.annotations['comment'][9:]

def get_myc_genomes():
  root = config.root
  gdir = os.path.join(root, 'genomes/Bacteria')

  gbks = []
  for root, dirs, files in os.walk(gdir):
    for f in files:
      if 'gbk' in f and 'Myco' in root:
        gbks.append(SeqIO.parse(os.path.join(root, f),'genbank').next())
        print 'extracting ' + f
  return gbks

      
  

def genomic_16s(nwk, ali, genomes):
  arch_list = []
  for g in genomes:
    org_name = g.annotations['organism'].replace(' ', '.')
    arch_list = nwk.root.clades[0].find_any(dict(name = '.*{0}.*'.format(org_name)))
    raise Exception()

def halo_16s(clade, ali):
  arch_list = clade.get_terminals()
  gbids_nwk = [clade_gbid(elt) for elt in arch_list]
  gbids_16s = [elt.annotations['comment'][9:] for elt in ali]

  ali_idxs = [gbids_16s.index(g) for g in gbids_nwk]
  #ali_recs = [ali[i] for i in ali_idxs]
  
  #slice to get a subalignment for species of interest.
  sub_a = None
  for a in ali_idxs:
    if sub_a == None:
      sub_a = ali[a:a+1]
    else:
      sub_a.append(ali[a])

  #and remove gaps:
  ungapped = 
