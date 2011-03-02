import compbio.projects.cbdb.adapters as adapters
from sqlalchemy import Column, Integer, String, Unicode, ForeignKey
import compbio.config as config

def get_maps():
  return dict(GBAcc = {'id': 2,
                      'accession':0,
                      'version':1})
                      
              
def get_tables():
  return [dict(name = 'GBAcc',
               attrs = {'id':Column(Integer,primary_key = True),
                       'accession':Column(String, index = True),
                       'version':Column(Integer)})]

def file_default():
  return config.dataPath('::unseen_data/gb_acclist.genbank')


def fill_db():
  filename = file_default()
  # a bit of a hack... there is only one table so it is just
  # the input filename
  table_defs = get_tables()
  column_maps =get_maps()
  table_names = map(lambda x: x['name'],table_defs)
  table_files = [filename]

  col_sep = ','
  rec_sep = '\n'

  adapters.fileMap2DB('gb_accid',
                      table_files,
                      table_names,
                      table_defs,
                      column_maps,
                      col_sep,
                      rec_sep,
                      reset = True)
                      
                      
