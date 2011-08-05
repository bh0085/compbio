import compbio.config as config
import re, itertools as it
import networkx as nx
import cb.p.network.io as nio
import cb.utils.memo as mem

def gene_symbol(fbid):
  gmap = open(config.dataPath('flybase/gene_map.tsv'))
  for l in gmap.xreadlines():
    if fbid + '\t' in l:
      return l.split('\t')[0]
  return ''
  
def gene_summary(fbid):
  gsums = open(config.dataPath('flybase/gene_summaries.tsv'))
  for l in gsums.xreadlines():
    if fbid in l[:len(fbid)]:
      return l
  return ''

def gene_biology(fbid):
  summ = gene_summary(fbid)
  try:
    bio = re.compile( ' ([^\.]*biolog[^\.]*\.)').search(summ).group(1)
  except Exception, e:
    return ''

  bio = re.sub(re.compile('<.*>'), '', bio)
  return bio

def gene_terms(fbid):
  bio = gene_biology(fbid)

  if 'process:' in bio:
    sub_str = bio.split('process:')[1]
    terms =  [l.strip() for l in sub_str.split(';')]
  elif 'under:' in bio:
    sub_str = bio.split('under:')[1]
    terms =  [l.strip() for l in sub_str.split(';')]
  elif bio == '':
    return set([])
  elif bio == 'The biological processes in which it is involved are not known.':
    return set([])
  else : 
    raise Exception()
  return set(terms)

def jacc(set1,set2):
  if len(set1) == 0 or len(set2) == 0: return 0.
  ''' returns the fraction of elements in set1 and set2 in both set1 and set2'''
  return float(len(set1.intersection(set2))) / len(set1.union(set1))
def spec(set1,set2):
  '''returns fraction of elements in set1 that are in set2'''
  return float(len(set1.intersection(set2))) / len(set1)
def sens(set1,set2):
  '''returns the fraction of elements of set2 covered by set1'''
  return float(len(set1.interesection(set2))) / len(set2)
def term_network_overlay(net_graph):
  simfunc = jacc 
  terms = dict([(gname, gene_terms(gname)) for gname in net_graph.nodes()])
  eints = dict([((e0,e1), simfunc(terms[e0],terms[e1])) 
                for e0,e1 in net_graph.edges()])
  
  return terms, eints


def term_groups(name = 'bdtnp', nterms = -1 ,**kwargs):
  '''
kwargs:
  nterms:  defaults to -1
'''
  def set_term_groups(**kwargs):
   nterms = kwargs.get('nterms') 
   if name == 'bdtnp': gene_list = nio.getBDTNP().keys()
   elif name == 'kn':    gene_list = graphs['kn'].nodes()
     
   #GET ALL CONTROLLED VOCAB TERMS APPLYING TO A GIVEN GENE LIST
   terms = [(gname,gt)  for gname in gene_list for gt in gene_terms(gname)
           ]
   all_terms = set([t[1] for t in terms])
   term_groups_tmp =[(k, list(g)) for k, g in 
                     it.groupby(
       sorted(terms, key = lambda x: x[1]),
       key = lambda x: x[1])
                     ]
   
   #SORT THE TERM GROUPS BY GENE COUNT AND ONLY TAKE TOP N
   if nterms == -1: nterms = len(term_groups_tmp)
   term_groups = sorted(term_groups_tmp, 
                        key = lambda x: len(x[1]))[::-1][:nterms]
   return term_groups
  return mem.getOrSet(set_term_groups, 
                      **mem.rc(kwargs,
                               on_fail = 'compute',
                               register = '{0}_{1}'.format(name, nterms),
                               nterms = nterms))

def term_network(name = 'bdtnp', nterms = -1 , **kwargs):
  '''
kwargs:
  nterms:  defaults to -1
'''
  
  def set_term_network( **kwargs):
    nterms = kwargs.get('nterms')
    name = kwargs.get('name')
    if name == 'bdtnp': gene_list = nio.getBDTNP().keys()
    elif name == 'kn':    gene_list = graphs['kn'].nodes()
    grps = term_groups(**mem.sr(kwargs,
                                name = name))
    network = nx.Graph()
    network.add_nodes_from(gene_list)
  
    
    for g in grps:
      edgelist =  [[g1[0],g2[0]]
                   for g1 in g[1] for g2 in g[1]]
      network.add_edges_from( edgelist  )
    return network
  return mem.getOrSet(set_term_network,
                      **mem.rc(kwargs,
                               on_fail = 'compute',
                               register = '{0}_{1}'.format(name,nterms),
                               nterms = nterms,
                               name = name))



