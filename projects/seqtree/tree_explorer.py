import compbio.projects.seqtree as sqt
import compbio.projects.cbdb as cbdb
from compbio.utils import pbar
from compbio.utils import gbid_lookup

import itertools as it

def terminal_maps(tree, reset = False):
  gba_dbi, tax_dbi, tgb_dbi = sqt.gbaDBI(), sqt.taxDBI(), sqt.tgbDBI()
  tmaps = {}
  for t in tree.get_terminals():
    acc = sqt.clade_gbacc(t)
    gbid = gba_dbi.Session.query(gba_dbi.GBAcc).filter_by(accession = acc).all()
    if not gbid: 
      print '{0:20} for nodename: {1}'.format('failure', t.name)
      tmaps[t.name] = None
      continue
    print '{0:20} for nodename: {1}'.format('sxs', t.name)

    node = ncbiNodeForGBID(gbid[0])
    sciname = [x for x in node.names if x.name_class == 'scientific name'][0].name_txt
    tmaps[t.name] = node.id

  return tmaps

def taxnode_accessions(nodes, 
                       querydb =  None):
  if querydb == None: querydb =cbdb.getName('group1.stk')
  print querydb
  tgb_dbi = cbdb.getName('tax_gbs')
  gba_dbi = cbdb.getName('gb_acc_idjoin')

  node_accessions = []
  bar = pbar.simple(len(nodes))
  for i in range(len(nodes)):
    n = nodes[i]
    bar.update(i)
    gbid_matches = [x.gbid 
                    for x in tgb_dbi.Session.query(tgb_dbi.TaxGBJoin).\
                      filter_by(taxid = n.id).all()]
    acc_matches = [[x.accession 
                   for x in gba_dbi.Session.query(gba_dbi.GBAcc).\
                     filter_by(id = gbid).all()]
                   for gbid in gbid_matches]
    node_accessions.append(acc_matches)
  bar.finish()

  ali_matches = [[querydb.Session.query(querydb.Sequence).\
                    filter_by(gb_accession = acc) for acc in node_ids]
                 for node_ids in node_accessions] 
  ali_matches = [list(it.chain(*node_ids)) for node_ids in node_accessions] 
  return ali_matches                
                    

def clade_accessions_forall(nodes):
  dblist = [cbdb.getName('hammerhead.stk'),
            cbdb.getName('group1.stk'),
            cbdb.getName('group2.stk')]

  alis=[ taxnode_accessions(nodes, db) for db in dblist ]
  return alis 
def draw_accessions():
  pass
    
def ncbiNodeForGBID(gbid):
  gba_dbi, tax_dbi, tgb_dbi = sqt.gbaDBI(), sqt.taxDBI(), sqt.tgbDBI()
  taxid =  tgb_dbi.Session.query(tgb_dbi.TaxGBJoin).\
      filter_by(gbid = gbid.id).one().taxid
  node = tax_dbi.Session.query(tax_dbi.Node).filter_by(id = taxid).one()    

  return node
