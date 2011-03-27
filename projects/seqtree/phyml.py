'''
Run phyml.

So far, only has one function: 
  phyml.tree(alignment, [runid])

Which returns a max likelihood tree given an alignment specified as a biopython Align.
'''

import compbio.config as config
import os
import Bio.Phylo as phylo
import Bio.AlignIO as aio
import subprocess

def tree(alignment,
         run_id = 'T%05i' % (0,)):
  old_cwd = os.getcwd()
  new_wd = config.dataPath('phyml')
  if not os.path.isdir(new_wd): os.mkdir(new_wd)
  os.chdir(new_wd)

  infilepath = 'infile{0}'.format(run_id)
  infile = open(infilepath,'w')
  aio.write(alignment, infile, 'phylip')
  infile.close()

  command = 'phyml -i {0} '.format(infilepath)
  print command
  subprocess.call(command,
                  shell = True)
  treefilepath = infilepath + '_phyml_tree.txt'
  treefile = open(treefilepath)
  tree =phylo.read(treefile, 'newick')
  treefile.close()
  os.chdir(old_cwd)
  return tree
  
