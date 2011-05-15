#!/usr/bin/env python
'''
paml.py:

Compute a maximum likelihood tree of ancestor sequences
given a preexisting tree in the biopython Newick "Tree"
format and a sequence alignment in the biopython "Align"
format. Optionally, give a run_id to save all files with
unique suffixes so that the algorithm may be safely run 
by multiple processes at once.

'''
import Bio
import Bio.Phylo as phylo
import Bio.AlignIO as aio
import compbio.config as config
import re, itertools as it, subprocess, os

def run_paml(tree_in,ali_in, run_id= 'T%05i' % (0,),
             verbose = False):
  '''
  Given an input tree in the form of a Biopython tree
  with branch lengths and names, write to a file and 
  run paml's baseml to generate a maximum likelihood 
  ancestry in the path data/paml/rst  '''

  paml_d = config.dataPath('paml')
  run_d = config.dataPath(os.path.join(paml_d , 'run_{0}'.format(run_id)))
  if not os.path.isdir(paml_d): os.mkdir(paml_d)
  if not os.path.isdir(run_d): os.mkdir(run_d)
  old_cwd = os.getcwd()
  os.chdir(run_d)

  outfilepath = 'paml_tree_{0}.paml'.format(run_id)
  
  treefilepath = 'paml_tree_{0}.newick'.format(run_id)
  treefile = open(treefilepath,'w')
  phylo.write(tree_in,treefile,'newick', plain = True)
  treefile.close()

  alifilepath ='paml_tree_{0}.phylip' .format(run_id)
  alifile = open(alifilepath, 'w')
  aio.write(ali_in, alifile, 'phylip')
  alifile.close()

  ctlfilepath= 'baseml_{0}.ctl'.format(run_id)
  ctlfile = open(ctlfilepath,'w')
  ctlfile.write(make_baseml(treefilepath,
                            alifilepath,
                            outfilepath,
                            ancestors = 1))
  ctlfile.close()

  command = 'baseml {0} '.format(ctlfilepath)

  #fix a damned paml bug.
  sed_command = "sed -i -e '1 s/$/\ \ I/' {0}"\
      .format(alifilepath)
  
  sprc = subprocess.Popen(sed_command, stdout = subprocess.PIPE, shell = True)
  comms = sprc.communicate()
  pprc = subprocess.Popen(    command, stdout = subprocess.PIPE, shell = True)
  comms = pprc.communicate()
  if verbose:
    print comms[0]

  os.chdir(old_cwd)
  rstfile = os.path.join(run_d,'rst')
  return rstfile
  
  
def make_baseml(treefile,
                alifile,
                outfile,
                ancestors = 1):
  '''Make a baseml.ctl file to use when running paml's baseml'''
  text = '''
seqfile = {alifile}
treefile = {treefile}
outfile = {outfile}    
RateAncestor = {ancestors}
verbose = 1
'''.format(alifile = alifile,
          outfile = outfile,
          treefile = treefile, 
          ancestors = ancestors)
  return text


def rst_parser(rstfile):
  #data = rstfile.read()
  rst = open(rstfile)
  #data = rst.read()
  

  
  while 1:
    if 'Ancestral reconstruction' in rst.readline(): break

  rst.readline()
  tree_actual_text =rst.readline()
  

  rst.readline()
  tree_symbolic_text =rst.readline()
  rst.readline()
  tree_edges_text = rst.readline()

  tree_terminal_nums = list([m for m in re.findall\
                              (re.compile('(\d+)'),\
                                 tree_symbolic_text)])
  tree_internal_nums = list(set([m for m in re.findall\
                                   (re.compile('(\d+)'),\
                                      tree_edges_text) 
                                 if not m in tree_terminal_nums]))


  rst.readline()
  rst.readline()
  tree_nodelabeled_txt =rst.readline()
  
  while 1:
    l = rst.readline()
    if 'site' in l and 'Freq' in l: break
  rst.readline()
  
  marginal_re = re.compile(\
    '\s*(?P<site>\d*)\s*(?P<Freq>\d*)\s*(?P<data0>\S*):\s*(?P<data>.*)')
  marginal_data_re = re.compile('([AGCT])\(([\d\.]*)\)')
  marginal = []
  while 1:
    l = rst.readline()
    init_parsing = marginal_re.search(l)
    if not init_parsing: break

    gdict = init_parsing.groupdict()
    gdict['data'] =[ x for x in 
                     marginal_data_re.findall(gdict['data'])]
    marginal.append(gdict)

  
  m_nums = [str(num) for num in sorted([int(n) for n in tree_internal_nums])]
  m_names= ['HN{0}'.format(item) for idx, item in enumerate(m_nums)]
  m_seqs =[''.join([m['data'][idx][0] for m in marginal])
           for idx in range(len(m_nums))]

  m_probs=[[float(m['data'][idx][1]) for m in marginal]
           for idx in range(len(m_nums))]
  
  all_edges = [m.groups() 
               for m in  re.compile('([\d]+)\.\.([\d]+)').\
                 finditer(tree_edges_text)]

               
  #it is important to note here that the terminal numbers are
  #listed in the same order as the terminal names on two seperate
  #lines of the file... 
  term_nums = tree_terminal_nums
  term_names = re.compile('(N\S*):').findall(tree_actual_text)
  if len(term_names) != len(term_nums):
    print 'OOPS... COULDNT PARSE TERM NAMES FROM RST.'
    raise Exception('Need to name nodes as N#*')
  
  m_clades = []
  for idx, m in enumerate(m_nums):
    node = phylo.Newick.Clade()
    node.name = m_names[idx]
    node_num = m_nums[idx]
    node.child_nums = [c for p, c in all_edges if p == node_num]
    node.m  = dict(seq = Bio.SeqRecord.SeqRecord(m_seqs[idx],
                                                 name = node.name),
                   probs = m_probs[idx])
    node.seq_probs = m_probs
    
    m_clades.append(node)
  
  t_clades = []
  for idx, t in enumerate(term_nums):
    node = phylo.Newick.Clade()
    node.name = term_names[idx]
    t_clades.append(node)

  all_children = set(it.chain(*[cl.child_nums for cl in m_clades]))
  m_root, =[ m_clades[m_nums.index(num)]
             for num in set(m_nums).difference(all_children)]
  
  for node in m_clades:
    children = [
      (m_clades + t_clades)[(m_nums+ term_nums).index(child)] 
      for child in node.child_nums]
    for c in children:
      node.clades.append(c)
       
    
  
  tree_out = phylo.Newick.Tree(m_root)
  #generate a tree from the bottom up:
  trees = {}
  terms = t_clades
      
       
  weighted_tree_in = phylo.NewickIO.Parser(tree_actual_text).\
      parse().next()
  terms_in = weighted_tree_in.get_terminals()
  for t0 in terms_in:
    tf, = [c for c in t_clades if c.name == t0.name]
    lens = [c.branch_length for c in weighted_tree_in.get_path(t0)]
    for idx, c in enumerate(tree_out.get_path(tf)): 
      c.branch_length = lens[idx]

  return tree_out
