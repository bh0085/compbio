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
  bio = re.compile( ' ([^\.]*biolog[^\.]*\.)').search(summ).group(1)
  bio = re.sub(re.compile('<.*>'), '', bio)
  return bio
