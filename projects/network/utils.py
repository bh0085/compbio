import compbio.config as config
import re
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
def term_network(net_graph):
  simfunc = jacc 
  terms = dict([(gname, gene_terms(gname)) for gname in net_graph.nodes()])
  eints = dict([((e0,e1), simfunc(terms[e0],terms[e1])) 
                for e0,e1 in net_graph.edges()])
  
  return terms, eints
