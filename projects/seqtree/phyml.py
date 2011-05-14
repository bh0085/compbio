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
         run_id = 'T%05i' % (0,),
         bionj = False):

  old_cwd = os.getcwd()
  new_wd = config.dataPath('phyml')
  if not os.path.isdir(new_wd): os.mkdir(new_wd)
  os.chdir(new_wd)

  infilepath = 'infile{0}'.format(run_id)
  infile = open(infilepath,'w')
  aio.write(alignment, infile, 'phylip')
  infile.close()


  command = 'phyml --quiet -i {0} -o {1} '.format(infilepath, 'n' if bionj else 'tlr' )
  print command
  subprocess.call(command,
                  shell = True,
                  stdout = subprocess.PIPE)
  treefilepath = infilepath + '_phyml_tree.txt'
  treefile = open(treefilepath)
  tree =phylo.read(treefile, 'newick')
  treefile.close()
  os.chdir(old_cwd)
  return tree
  
