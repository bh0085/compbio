import re
import inspect
import subprocess
import compbio.utils.path_mgr as pm
from compbio.utils.path_mgr import *
from Bio import SeqIO


def testIO():
  fname = 'genomes/Bacteria/Mycoplasma_mycoides_SC_PG1_uid58031/NC_005364.fna'
  fin = open(fname)

  pushd()
  read = SeqIO.parse(fin, 'fasta')
  pushd('temp')
  
  fout =open('temp.clw','w')
  SeqIO.write(read, fout, 'clustal')
  fin.close()
  fout.close()
  popd()
  popd()

    
def load(filename = 'data/rRNA'):


  fname = 'genomes/Bacteria/Mycoplasma_mycoides_SC_PG1_uid58031/NC_005364.fna'
  fin = open(fname)

  pushd()
  read = SeqIO.parse(fin, 'fasta')
  pushd('temp')  
  cfile = 'temp.clw'
  fout =open(cfile,'w')
  SeqIO.write(read, fout, 'clustal')
  fin.close()
  fout.close()
  
  f = open(cfile)
  seq = SeqIO.parse(cfile,'clustal')
  seq_l = list(seq)
  seq_str = seq_l[0].seq.tostring()
  
  l = len(seq_str)
  ofs = 200
  n = l / ofs - 1

  for i in range(n):
    prc = subprocess.call('RNAz --window={0}-{1} --output {2} temp.clw '.format( i * ofs, ( i + 1 * ofs), 'z.rnaz'),
                          shell = True)
    z = open('z.rnaz').read()
    print z
    break

  pm.popd()
  pm.popd()
  
