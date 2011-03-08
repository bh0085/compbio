import compbio.projects.seqtree as sqt
from compbio.projects.cbdb import gb_alignment as gb_ali
import compbio.projects.cbdb as cbdb
from compbio.utils import pbar
from compbio.utils import gbid_lookup

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

def taxnode_accessions(nodes):
  
  ali_dbi = cbdb.getName('hammerhead.stk', gb_ali.get_tables())
  tgb_dbi = sqt.tgbDBI()
  gba_dbi = sqt.gbaDBI()

  node_accessions = {}
  bar = pbar.simple(len(nodes))
  for i in range(len(nodes.keys())):
    n = nodes.values()[i]
    bar.update(i)

    gbid_matches = [x.gbid for x in tgb_dbi.
Session.query(tgb_dbi.TaxGBJoin).filter_by(taxid = n).all()]
    acc_matches = []
    for gbid in gbid_matches:
      gbacc_matches = ['{0}.{1}'.format(x.accession,x.version) for x in gba_dbi.Session.query(gba_dbi.GBAcc).filter_by(id = gbid).all()]
      
      acc_matches.extend(gbacc_matches)
      
    node_accessions[nodes.keys()[i]] =acc_matches


  bar.finish()
  print 
  print 'finding accession matches in rfam family'
  matches ={}
  ct = 0 
  bar = pbar.simple(len(nodes))
  for k,v in node_accessions.iteritems():
    
    bar.update(ct)
    ct+= 1
    matches[k]=[]

    for acc in v:
      g2match = [x.id for  x in ali_dbi.Session.query(ali_dbi.Sequence).filter_by(gb_accession=acc).all()]
      if g2match: matches[k].extend(g2match)
  bar.finish()
  return matches, node_accessions
                    
                    
    
def ncbiNodeForGBID(gbid):
  gba_dbi, tax_dbi, tgb_dbi = sqt.gbaDBI(), sqt.taxDBI(), sqt.tgbDBI()
  taxid =  tgb_dbi.Session.query(tgb_dbi.TaxGBJoin).filter_by(gbid = gbid.id).one().taxid
  node = tax_dbi.Session.query(tax_dbi.Node).filter_by(id = taxid).one()    

  return node
