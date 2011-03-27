#!/usr/bin/env python
'''
The btol class contains an instances of the microbial tree of life annotated
with nodes in the ncbi taxonomy whenever nodes are available for the genbank
accessions with which terminals are marked.

BTOL contains one batch processing script at the moment

BTOL also contains methods with which to map non-terminal nodes to phyla when
possible, to fetch sequences from mapped alignment databases and put them
at appropriate nodes and to fetch genomes as well.

--e.g: set_ranks(rank = 'phylum')
       set_seqs(dbname = 'hammerhead')
       set_genomes()

Using the metric on the tree, it is able to compute distances between sequence

--e.g: get_dist(node1, node2)

'''
import inspect
import itertools as it, os
import compbio.utils.pbar as pbar

import Bio.Phylo.Newick as nt
from Bio.Seq import *
import Bio.Phylo as Phylo

from compbio.projects.seqtree import muscle, phyml, paml
import compbio.projects.cbdb as cbdb
from compbio.projects.cbdb.helpers import ncbi, ali 
import compbio.projects.seqtree as sqt
import compbio.utils.memo as mem
import compbio.utils.bsub as bsub
from numpy import *
import matplotlib.pyplot as plt


def run(**kwargs):
  BT = getBTOL(**mem.sr(kwargs))
  seqnodes = BT.investigatePhylum(**kwargs)
  recs, seqelts, seqtuples = seq_recs(seqnodes)
  align = align_seqnodes(recs)
  tree = phyml.tree(align)
  rstfile= paml.run_paml(tree, align)
  anc_tree = paml.rst_parser(rstfile)

  anc_alignment = [SeqRecord(elt.m['seq'], 
                             id = None,
                             name = elt.name,
                             annotations = {'scores':elt.m['probs']})
                   for elt in anc_tree.get_nonterminals()]
  

  return (tree, anc_tree), (align, anc_alignment)


def seq_recs(seqnodes, **kwargs):
  seqtuples = sorted([
      (sn.seqs[k].replace('-',''),sn.id, sn.seqids[k], k)
      for sn in seqnodes if len(sn.seqs) > 0  
      for k  in sn.seqs ], key = lambda x: x[0])
  
  seqletters, seqelts = zip(*[(k,list(g)) 
                              for k, g in it.groupby(seqtuples, key = lambda x: x[0])])
  recs = [ (SeqRecord(Seq(elt), id = str(seqelts[i][0][3]), name = seqelts[i][0][1]))
            for i,elt in enumerate(seqletters) ]
  
  return recs,  seqelts, seqtuples
def align_seqnodes(recs,**kwargs):
  align = muscle.align(recs[0:5])
  return align


def getBTOL(**kwargs):
  def setBTOL(**kwargs):
    B = BTOL(**mem.sr(kwargs))
    if not B.treeInitialized():
      print 'Underlying tree structure apparently uninitialized: initializing\n...'
      B.initTree()
      print '...\nDone\nSaving\n...'
      B.saveTree()
      print '...\nDone'
    return B
  return mem.getOrSet(setBTOL, **mem.rc(kwargs, register = 'BTOL'))

class BTOL():
  '''BTOL Class - see module doc: btol?'''
  def __init__(self,
               **kwargs):
    self.t  = sqt.init(**mem.sr(kwargs))



  def treeInitialized(self):
    if not 'phylum' in it.chain(*[ leaf.m.keys() for leaf in self.t.get_terminals()[0:10]]):
      return False
    else: return True
  def initTree(self):
    self.makeRank('phylum')
    self.makeRank('family')
      
  def saveTree(self):
    sqt.init(update = self.t)
  def makeRank(self, rank = 'phylum', subtree = None): 
    #Get the subtree and db connections to build meta for
    tree = subtree if subtree != None else self.t     
    dbi = cbdb.getName('taxdmp')

    print 'Fetching taxonomic nodes from the db'
    #Get the terminal nodes and corresponding ncbi taxa
    terms = [t for t in tree.get_terminals() if t.m.has_key('taxid')]
    nodes = [dbi.S.q(dbi.Node).filter_by(id = t.m['taxid']).scalar() 
             for t in terms]

    #endpoints for parental iteratiion
    taxa = ncbi.get_rank(rank)
    root = ncbi.get_root()


    print 'Computing terminal node mappings for taxon: {0}'.format(rank)
    bar = pbar.simple(len(nodes)); bar.start()

    node_taxa = list(nodes)
    get_p_iter = lambda: \
        node_taxa[idx] == None and True \
        or node_taxa[idx] in taxa and True \
        or node_taxa[idx] == root and True \
        or node_taxa.__setitem__(idx,node_taxa[idx].parent) \
        or node_taxa[idx]

    for idx, v in enumerate(node_taxa):
      bar.update(idx);
      par = list(iter(get_p_iter, True))[-1] if v else None
      terms[idx].m[rank] = par.id if par in taxa else None
    bar.finish()
    print 'Done!'

  
                       
  def leafNodes(self,**kwargs):
    def setLeafNodes(**kwargs):
      all_leaves = self.t.get_terminals()
      dbi = cbdb.getName('taxdmp')
      all_nodes = [ ncbi.get_node(l.m['taxid'],dbi) 
                    if 'taxid' in l.m.keys() else None for l in all_leaves]
      return all_nodes
    nodes = mem.getOrSet(setLeafNodes, 
                         **mem.rc(kwargs,
                                  hardcopy = False,
                                  on_fail = 'compute',
                                  register = 'leaf_nodes'))
    return nodes

  def getTaxon(self,rank = rank,
               **kwargs):
    def setTaxon(BTInstance = None, rank = None, **kwargs):
      assert rank; assert BTInstance
      leafnodes = BTInstance.leafNodes(**mem.sr(kwargs))
      leaf_families = [ncbi.get_taxon(node, rank=rank) 
                       if node else None for node in leafnodes]
      return leaf_families
    return mem.getOrSet(setTaxon,
                        **mem.sr(kwargs, 
                                 rank = rank, 
                                 BTInstance = self,
                                 on_fail = 'compute',
                                 hardcopy = False,
                                 register = rank))

  def ali_in_tree(self,aliname = 'group2.stk',
                  rank = 'genus',
                  **kwargs):
    all_seqs = ali.get_seqs(aliname)
    alinodes = ali.get_taxnodes(aliname)

    aliranks = [t.rank if t else None for t in alinodes]
    all_leaves = self.t.get_terminals()
    leafnodes = self.leafNodes(reset = mod(reset, 2))
    leafranks =[n.rank if n else None for n in leafnodes]
 
    ali_families = ali.get_taxon_forall(rank = rank,aliname = aliname,
                                        **mem.sr(kwargs))
    leaf_families=  self.getTaxon(rank = rank, **mem.sr(kwargs))

    aset = set(ali_families)
    lset = set(leaf_families)

    a_domains =[(node, ncbi.get_taxon(node,'superkingdom')) for node in aset]
    l_domains =[(node, ncbi.get_taxon(node,'superkingdom')) for node in lset]
    
    bac_domain = [x[1] for x in l_domains if ncbi.sciname(x[1])== 'Bacteria'][0]

    l_bacs = set((l[0] for l in l_domains if l[1] == bac_domain))
    a_bacs = set((a[0] for a in a_domains if a[1] == bac_domain))

    leaf_bacteria = [leaf  if leaf in l_bacs else None for leaf in leaf_families]
    ali_bacteria =  [a  if a  in a_bacs else None for a in ali_families]
  
    return leaf_bacteria, ali_bacteria, leafnodes, alinodes

  def show_rank(self,
                tree_n_rank, ali_n_rank,
                tree_n, ali_n
                ):


    tree_leaves =self.t.get_terminals()

    rank_unq = set(tree_n_rank).union(set(ali_n_rank))
    supertaxa = {}
    for u in rank_unq:
      if not u: continue
      
      tree_matches = nonzero(equal(tree_n_rank, u))[0]
      tree_leaf_subset = [tree_leaves[i] for i in tree_matches]
      tree_n_subset = [tree_n[i] for i in tree_matches]

      ali_matches = nonzero(equal(ali_n_rank, u))[0]
      ali_n_subset = [ali_n[i] for i in ali_matches]

      supertaxa[ncbi.sciname(u)] = (tree_n_subset, ali_n_subset)

    f = plt.figure(0)  
    f.clear()
    xax = arange(len(supertaxa))
    ax = f.add_subplot('111',xticks = [])
    
    plots = [ax.plot(xax,[log(len(e[0])+1) for e in supertaxa.values()],
                     linewidth = 2, label='tree')[0],
             ax.plot(xax,[log(len(e[1])+1) for e in supertaxa.values()],
                     linewidth = 2, label='ali')[0]]
    ax.legend(plots, [x.get_label() for x in plots])
    for i,x in enumerate(xax):
      ax.text(x,0,supertaxa.keys()[i], 
              va = 'top', ha = 'left',
              rotation = -15.)

  
  def investigatePhylum(self, 
                        aliname = 'group2.stk',
                        p_node = None, **kwargs):

    if not p_node: p_node = ncbi.taxon_with_name('phylum', 'Thermotogae')
    ali_seqs = ali.get_seqs(aliname, **mem.sr(kwargs))
    ali_nodes = array(ali.get_taxnodes(aliname, **mem.sr(kwargs)))
    ali_phyla = array(ali.get_taxon_forall(aliname,**mem.sr(kwargs, rank = 'phylum')))
    ali_inds = nonzero(equal(ali_phyla, p_node))[0]
    
    leaf_terminals = self.t.get_terminals()
    leaf_nodes = array(self.leafNodes(**mem.sr(kwargs)))
    leaf_phyla = array(self.getTaxon('phylum', **mem.sr(kwargs)))
    leaf_inds = nonzero(equal(leaf_phyla, p_node))[0]
        
    ap_sub = ali_phyla[ali_inds]
    lp_sub = leaf_phyla[leaf_inds]

    ag_sub = array(ali.get_taxon_forsome(ali_nodes[ali_inds],'genus','thermo',
                                         **mem.sr(kwargs)))
    lg_sub = array(self.getTaxon('genus', **mem.sr(kwargs)))[leaf_inds]

    as_sub = array(ali.get_taxon_forsome(ali_nodes[ali_inds], 'species', 'thermo'))
    ls_sub = array(self.getTaxon('species',**mem.sr(kwargs)))[leaf_inds]

    db16 = cbdb.getName('16s')
    a_16s= [ db16.S.q(db16.Sequence).
             filter_by(source_taxon = n.id).all() for n in ali_nodes[ali_inds]]
    l_16s= [ db16.S.q(db16.Sequence).
             filter_by(source_taxon = n.id).all() for n in leaf_nodes[leaf_inds]]

    all_lens = dict([ (k, [len(list( e)) for e in seqlist] )  
                      for seqlist,k in [[a_16s,'a_16s'],[l_16s,'l_16s']]])
    

    leaf_sns = [ SeqNode(lg_sub[i],ls_sub[i] , leaf_nodes[idx], 
                         [(x.sequence,x.gb_id) 
                          for x in l_16s[i]],
                         src = leaf_terminals[idx],
                         node_id =  'btol:default:{0}'.format(leaf_terminals[idx].m['id']))
                 for i, idx in enumerate(leaf_inds)]
    ali_sns = [ SeqNode(ag_sub[i],as_sub[i] , ali_nodes[idx],
                        [( x.sequence,x.gb_id ) 
                         for x in a_16s[i]],
                        src = ali_seqs[idx], 
                        node_id = 'ali:{0}:{1}'.format(aliname,ali_seqs[idx].id))
                for i, idx in enumerate(ali_inds)]

    return list(it.chain(leaf_sns, ali_sns))


sn_ct = 0
class SeqNode(object):
  def __init__(self, 
               family_node,
               species_node,
               term_node,
               seqs,
               src = '',
               node_id = '::'):
    self.family_node = family_node
    self.term_node = term_node
    self.species_node = species_node
    self.src = src
    self.id = node_id

    global sn_ct
    self.name = 'N%03i' % (sn_ct,)
    sn_ct += 1

    self.seqs = dict([(self.name+'_%02i'%idx , seq[0]) for idx, seq in enumerate(seqs)])
    self.seqids = dict([(self.name+'_%02i'%idx , seq[1]) for idx, seq in enumerate(seqs)])


  def __repr__(self):
    return '{0},{1},{2} ({3})'.format(ncbi.sciname(self.family_node),
                                      ncbi.sciname(self.species_node),
                                      ncbi.sciname(self.term_node),
                                      self.term_node.rank)

    

def run_anc(input_dict,run_id = None):
  assert run_id
  rank_name = input_dict['rank_name']
  taxon_id  = input_dict['taxid']
  aliname = input_dict['aliname']
  
  BT = getBTOL()
  p_node = ncbi.get_node(taxon_id)
  seqnodes = BT.investigatePhylum(p_node = p_node)
  recs, seqelts, seqtuples = seq_recs(seqnodes)
  align = align_seqnodes(recs)
  tree = phyml.tree(align, run_id = run_id)
  rstfile= paml.run_paml(tree, align, run_id = run_id)
  anc_tree = paml.rst_parser(rstfile)

  anc_alignment = [SeqRecord(elt.m['seq'], 
                             id = None,
                             name = elt.name,
                             annotations = {'scores':elt.m['probs']})
                   for elt in anc_tree.get_nonterminals()]
  out_dict = dict(anc_tree=anc_tree,
                  anc_align= anc_alignment,
                  term_tree = tree,
                  term_align = align)
  return out_dict
  
def make_anc_batches(aliname, rank_name = 'phylum', 
                     do_bsub = False,
                     run = True, nrun = 0, 
                     **kwargs):
  BT = getBTOL(**mem.sr(kwargs))
  tree_tax = BT.getTaxon(rank_name)
  ali_tax = ali.get_taxon_forall(aliname, rank_name)
  union_tax = set(tree_tax).intersection(set(ali_tax))
  union_taxids = [e.id for e in union_tax if e]

  batch_pdicts = [dict(taxid = taxid,
                       rank_name = rank_name,
                       aliname = aliname)
                  for taxid in union_taxids]


  cmds = []
  for idx,  d in enumerate(batch_pdicts):
    run_id=bsub.get_run_id(idx, prefix = rank_name)
    bsub.save_inp(d, run_id)
    cmds.append(bsub.cmd(os.path.abspath(inspect.stack()[0][1]),'run_anc',run_id, do_bsub = do_bsub, run_id = run_id))

  #if run:
    
  return cmds
               

def usage():
  print '''
usage: btol.py run_id

Run a batch process on BTOL with inputs stored in 
data/batch/inputs/{run_id}.inp in pickle serial.
'''
if __name__ == "__main__":
  if len(sys.argv) < 3: usage()
  run_id = sys.argv[2]
  run_func = globals()[sys.argv[1]]
  input_dict = bsub.load_inp(run_id)
  output_dict = run_func(input_dict, run_id)
  bsub.save_out( output_dict, run_id)

  exit(0)



