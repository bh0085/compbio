from compbio.utils import path_mgr as pm, memo as mem
import os
from Bio import AlignIO

def myco(  reset = False):'''
For now, leave this project on hold until matt tells me an easy way to line up all elements of the family on a family tree...
'''
  myd  = pm.getd()

  if reset:
    allfile = os.path.join(myd,'RF00029.full')
    seedfile = os.path.join(myd,'RF00029.seed')
    
    s_seqs = AlignIO.parse(seedfile,'stockholm').next() 
    a_seqs = AlignIO.parse(allfile,'stockholm').next()

    mem.write(register = 's_seqs', value = s_seqs)
    mem.write(register = 'a_seqs', value = a_seqs)
  else:
    s_seqs, sxs1 = mem.read(register = 's_seqs', hardcopy = False)
    a_seqs, sxs2 = mem.read(register = 'a_seqs', hardcopy = False )
    if sxs1 + sxs2 != 2:
      raise Exception()
    
  for s in s_seqs:
    print 'id: {0:100}'.format(s.id)

  for a in a_seqs[0:100]:
    print 'id: {0:100}'.format(a.id)

  for a in a_seqs:
    if 'SACH' in a.id.upper():
      print a.id

  
  
  
