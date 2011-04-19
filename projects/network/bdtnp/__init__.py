import parser 
import exp
import os
from compbio import config
all = ['parser', 'exp']
if not os.path.isdir(config.dataPath('bdtnp')):
  os.mkdir(config.dataPath('bdtnp'))
