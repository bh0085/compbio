#!/usr/bin/env python

import compbio.config as config
import os
import re

def sort_prefixes(volume_name = 'cb'):
  prefix_path = config.dataPath(config.dataURL('genbank/prefixes'))
  for p in os.listdir(prefix_path):
    f= os.path.join(prefix_path, p)
    fopen = open(f)
    lines = fopen.readlines()
    lsort = sorted(lines)
    fopen.close()
    fopen = open(f ,'w')
    fopen.writelines(lsort)
    fopen.close()
    print p
    
def prefix(query):
  return re.compile('[A-Z]*').search(query).group()
def search_sorted(prefix_name, query, volume_name = 'cb'):
  '''performs a binary search within a sorted file to find the
  genbank id for a given query'''

  prefix_file = os.path.join(\
    config.dataPath(\
      config.dataURL('genbank',
                     volume_name = volume_name)),
    'prefixes/'+prefix_name)
  fopen = open(prefix_file)
  size = os.path.getsize(prefix_file)
  start = 0
  stop = size
  
  hplast = 0
  while 1:
    halfpt = (start + stop) / 2    
    fopen.seek(halfpt)
    if halfpt == 0:
      line = fopen.readline()
    else:
      blank = fopen.readline()
      line = fopen.readline()

    c0 = line.split(',')[ 0]

    if c0 == query: return line.split(',')[2].strip()
    elif c0 < query: 
      start = halfpt
    else:
      stop = halfpt
    
    if halfpt == hplast: raise Exception('Query not for: %s'%query)
    hplast = halfpt
    if start == stop: raise Exception('Query not found: %s' % query)
      
    
def split_prefixes(volume_name = 'cb'):
  '''splits the massive genebank accession list up by prefixes.
takes a volume name as a parameter in case the accesion list is
stored in an atypical location'''

  path = config.dataPath(config.dataURL('genbank/gb_acclist.genbank', 
                                        volume_name =volume_name)
                         )
  path_home = os.path.dirname(path)
  fopen = open(path)
  prefixes = {}
  count = 0
  long_count = 0
  for l in fopen.xreadlines():
    if l[1].isdigit(): pend = 1
    elif l[2].isdigit(): pend = 2
    elif l[3].isdigit(): pend = 3
    elif l[4].isdigit(): pend = 4
    elif l[5].isdigit(): pend = 5
    elif l[6].isdigit(): pend = 6
    else: raise Exception()

    prefix = l[0:pend]
    
    if not prefixes.has_key(prefix):
      print 'getting pre'
      prefixes[prefix] =  open(os.path.join(path_home, 'prefixes/'+prefix), 'a')
    f = prefixes[prefix]
    f.write(l)

    count += 1
    if count > 100000:
      count = 0
      long_count += 1
      print prefix, l
      if long_count > 10:
        print prefixes
        long_count = 0 
        while prefixes:
          f = prefixes.pop(prefixes.keys()[0])
          f.close()


  for k, p in prefixes.iteritems():
    p.close()

  
if __name__ == '__main__':
  '''by default, parses the full genbank accession list into seperate files for each prefix'''
  split_prefixes()
  exit()
