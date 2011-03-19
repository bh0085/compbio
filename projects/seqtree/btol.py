import Bio.Phylo.Newick as nt
import compbio.projects.seqtree as sqt

class BTOL(nt.Tree):
  def __init__(self,
               reset = False):
    self.t  = sqt.init(reset = False)
    #raise Exception()

  def list_phyla(self):
    pass
