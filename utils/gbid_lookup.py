#!/usr/bin/env python

import compbio.config as config
import os

def split_prefixes(volume_name = 'cb'):
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
      prefixes[prefix] =  open(os.path.join(path_home, prefix), 'a')
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
  split_prefixes()
  exit()
