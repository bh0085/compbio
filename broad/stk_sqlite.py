import os
import compbio.config as config
for r, d, fs in os.walk(config.dataPath(config.dataURL('unseen_data'))):
  for f in fs:
    if '.stk' in f:
      print f
