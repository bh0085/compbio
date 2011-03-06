from sqlalchemy import Column, Integer, String, Unicode, DateTime, ForeignKey, UniqueConstraint
from sqlalchemy.orm import relation
from compbio.projects import cbdb
from compbio import config
import os
import numpy as np


def get_maps():
  return dict(Node = {'id': 0,
                      'parent_taxid':1,
                      'rank':2,
                      'embl_code':3,
                      'gencodeid':6,
                      'mit_gencodeid':8,
                      'comments':12},
              Name = {'nodeid':0,
                      'name_txt':1,
                      'unique_name':2,
                      'name_class':3},
              Gencode = {'id':0,
                         'abbreviation':1,
                         'name':2,
                         'cde':3,
                         'starts':4},
              Citation = {'id':0,
                          'key':1,
                          'pubmed_id':2,
                          'medline_id':3,
                          'url':4,
                          'text':5,
                          'taxid_list':6}
              )


def get_tables():
  return [dict(name = 'Node',
               attrs = {'id':Column(Integer,primary_key = True),
                        'parent_taxid': Column(Integer, ForeignKey('node.id')),
                        'parent':relation('Node',
                                          primaryjoin="Node.parent_taxid==Node.id",
                                          backref = 'children'),
                        'rank': Column(String),
                        'embl_code':Column(String),
                        'mit_gencodeid':Column(Integer,ForeignKey('gencode.id')),
                        'gencodeid':Column(Integer, ForeignKey('gencode.id')),
                        'gencode':relation('Gencode', 
                                           primaryjoin="Gencode.id==Node.gencodeid"
                                           ),
                        'mit_gencode':relation('Gencode', 
                                      primaryjoin="Gencode.id==Node.mit_gencodeid"
                                               ),
                        'comments':Column(String)
                       }),
          dict(name = 'Name',
               attrs = {'nodeid':Column(Integer, ForeignKey('node.id')),
                        'node':relation('Node',backref = 'names'),
                        'name_txt':Column(String),
                        'unique_name':Column(String),
                        'name_class':Column(String)
                       }),
          dict(name = 'Gencode',
               attrs = {'id':Column(Integer,primary_key=True),
                       'abbreviation':Column(String),
                       'name':Column(String),
                       'cde':Column(String),
                       'starts':Column(String)}),
          dict(name = 'Citation',
               attrs = {'id':Column(Integer, primary_key = True),
                        'key':Column(String),
                        'pubmed_id':Column(Integer),
                        'medline_id':Column(Integer),
                        'url':Column(String),
                        'text':Column(String),
                        'taxid_list':Column(String)})]
def fill_db( reset = False):
  dbi = cbdb.getName('taxdmp', tables = get_tables(), reset = np.mod(reset, 2))
  filepath = config.dataPath('ncbi/taxdmp')
  fsize = os.path.getsize(filepath)
  maps = get_maps()

  record_sep = '\t|\n'
  col_sep = '\t|\t'
  colfun = lambda x: unicode(x, errors = 'replace').replace(record_sep, '').split(col_sep)
  record_iterfun = lambda x: x.xreadlines()
  
  fill_tables = {'Gencode':'gencode.dmp',
                 'Node':'nodes.dmp',
                 'Name':'names.dmp',
                 'Citation':'citations.dmp'}

                 
  count = 0
  for k,v in fill_tables.iteritems():
    fopen = open(os.path.join(filepath, v))
    mapped_class = dbi.__dict__[k]
    mapped_columns = maps[k]
    l0 = ''
    for l in record_iterfun(fopen):
      count += 1
      l0+=l
      if l0[-3:] == record_sep :
        l = l0
        l0 = ''
      else: continue
      cols = colfun(l)
      cls = mapped_class(**dict(map(lambda (x,y): (x,cols[y]),
                                    mapped_columns.iteritems())))
      dbi.Session.merge(cls)
      if np.mod(count, 1000) == 0:
        dbi.Session.commit()
        print k, v, count, cols, '{0:4}%'.format(100 * float(fopen.tell()) / fsize)
    dbi.Session.commit()
  return

if __name__ == '__main__':
  print 'filling db'
  fill_db()
  exit(0)

                    
      
'''   
citations.dmp
delnodes.dmp
division.dmp
gencode.dmp
merged.dmp
names.dmp
nodes.dmp

nodes.dmp
---------

This file represents taxonomy nodes. The description for each node includes 
the following fields:

	tax_id					-- node id in GenBank taxonomy database
 	parent tax_id				-- parent node id in GenBank taxonomy database
 	rank					-- rank of this node (superkingdom, kingdom, ...) 
 	embl code				-- locus-name prefix; not unique
 	division id				-- see division.dmp file
 	inherited div flag  (1 or 0)		-- 1 if node inherits division from parent
 	genetic code id				-- see gencode.dmp file
 	inherited GC  flag  (1 or 0)		-- 1 if node inherits genetic code from parent
 	mitochondrial genetic code id		-- see gencode.dmp file
 	inherited MGC flag  (1 or 0)		-- 1 if node inherits mitochondrial gencode from parent
 	GenBank hidden flag (1 or 0)            -- 1 if name is suppressed in GenBank entry lineage
 	hidden subtree root flag (1 or 0)       -- 1 if this subtree has no sequence data yet
 	comments				-- free-text comments and citations

names.dmp
---------
Taxonomy names file has these fields:

	tax_id					-- the id of node associated with this name
	name_txt				-- name itself
	unique name				-- the unique variant of this name if name not unique
	name class				-- (synonym, common name, ...)

division.dmp
------------
Divisions file has these fields:
	division id				-- taxonomy database division id
	division cde				-- GenBank division code (three characters)
	division name				-- e.g. BCT, PLN, VRT, MAM, PRI...
	comments

gencode.dmp
-----------
Genetic codes file:

	genetic code id				-- GenBank genetic code id
	abbreviation				-- genetic code name abbreviation
	name					-- genetic code name
	cde					-- translation table for this genetic code
	starts					-- start codons for this genetic code

delnodes.dmp
------------
Deleted nodes (nodes that existed but were deleted) file field:

	tax_id					-- deleted node id

merged.dmp
----------
Merged nodes file fields:

	old_tax_id                              -- id of nodes which has been merged
	new_tax_id                              -- id of nodes which is result of merging

citations.dmp
-------------
Citations file fields:

	cit_id					-- the unique id of citation
	cit_key					-- citation key
	pubmed_id				-- unique id in PubMed database (0 if not in PubMed)
	medline_id				-- unique id in MedLine database (0 if not in MedLine)
	url					-- URL associated with citation
	text					-- any text (usually article name and authors)
						-- The following characters are escaped in this text by a backslash:
						-- newline (appear as "\n"),
						-- tab character ("\t"),
						-- double quotes ('\"'),
						-- backslash character ("\\").
	taxid_list				-- list of node ids separated by a single space
'''

