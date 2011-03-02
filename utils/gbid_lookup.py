import compbio.config as config
import os

def lookup(accession, 
           data_volume = 'cb', depth = 20):
  path = config.dataPath(config.dataURL('genbank/gb_acclist.genbank', 
                                        volume_name = data_volume)
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
    prefixes.get(prefix, open(os.path.join(path_home, prefix), 'w')).write(l)
    if len(prefixes.keys()) > 5:
      break
    count += 1
    if count > 1000:
      count = 0
      long_count += 1
      print prefix
      if long_count > 10:
        long_count = 0 
        while prefixes:
          f = prefixes.pop(prefixes.keys()[0])
          f.close()


  for k, p in prefixes.iteritems():
    p.close()

  
