from sqlalchemy import Column, Integer, String, Unicode, DateTime, ForeignKey, UniqueConstraint, MetaData, create_engine
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relation, scoped_session, sessionmaker
from compbio import config
import os
import numpy as np
import itertools as it
from sqlalchemy.exc import OperationalError
import psycopg2
from sqlalchemy.sql import text

#import gb_alignment, gb_genomes, ncbi_tax, ncbi_tax_gbjoin
__all__ = [os.path.splitext(f)[0] 
           for f in os.listdir(os.path.dirname(globals()['__file__'])) #
           if os.path.splitext(f)[-1] == '.py' and f[0:2] != '__']
#__all__ = ['ncbi_tax']

def tableDefs(name):
  return globals()[name].get_tables()

def makeWithNameAndTables(name,  tables,
                          postgres = False,
                          reset = 0, host = None):

  #Connect an engine to the database
  engine = connectDB(name, postgres = postgres, host = host)

  #And generate a metadata object for the tables provided
  base, table_list = metadataWithTables(tables)
  Session = scoped_session(sessionmaker(autocommit = False, 
                                        bind = engine))
  
  #Connect to the dbfile; if it does not exist then create one.
  if reset:
    #inp = raw_input('Really delete existing database for "{0}"? [y/n]'.format('hi'))
    #if inp != 'y': raise Exception('Insert cancelled - either permit deletion or rerun with "reset = 0"')
    print 'Deleting the current database for %s' % name
    
    print postgres
    if not postgres:
      dbtables = map(lambda x: x[0], Session.execute('''SELECT name FROM sqlite_master WHERE type='table';''').fetchall())
      for t in base.metadata.tables:
        if t in dbtables:
          Session.execute('drop table %s;' % t)
          Session.commit()
      base.metadata.create_all(bind = engine, checkfirst = False)


    else:
      dbtables =  map(lambda x: x[0], Session.execute('select tablename from pg_tables').fetchall())
      for t in base.metadata.tables:
        if t in dbtables:
          Session.execute('drop table %s cascade;' % t)
          Session.commit()
      base.metadata.create_all(bind = engine, checkfirst = False)

 

  #Create a dbi storing the session and table types.xs
  attrs = dict(Session = Session, 
               base = base,
               **table_list)
  dbi = type(name+'DBI', (object,), attrs)
  globals()[name] = dbi

#return an engine connected to the name
def connectDB(name, postgres =False, **kwargs):
  print postgres
  if postgres:
    try:
      engine = create_engine(config.postgresPath(name,**kwargs))
      cxn = engine.connect()
      cxn.close()
    except OperationalError, e:
      conn = create_engine(config.postgresDefault(**kwargs)).connect()
      conn.connection.set_isolation_level(psycopg2.extensions.ISOLATION_LEVEL_AUTOCOMMIT)
      conn.execute("""create database  """+  name+ ";")
      conn.close()
      engine = create_engine(config.postgresPath(name,**kwargs))
    return engine
  else:
    print config.sqlitePath(name)
    engine = create_engine(config.sqlitePath(name))
    return engine


def getName(name,tables = None, postgres = False, host = None, reset = 0):
  if reset or not name in globals().keys():
    makeWithNameAndTables(name,postgres = postgres, tables = tables, reset = np.mod(reset,2),host = host)
  return globals()[name]

def metadataWithTables(tables):
  meta = MetaData()
  tlist = {}

  base =  declarative_base(metadata = meta)
  for t in tables:
    name = t['name']
    attrs = t['attrs']
    tlist = dict(_initDeclarative(name, attrs, base), **tlist)

  return base, tlist

def _initDeclarative(name, attrs_in, base):
  attrs = dict(attrs_in)
  if not '__tablename__' in attrs.keys():
    attrs['__tablename__'] = str.lower(name)
  
    
  #check if the table has a primary key.
  if not len(np.nonzero(map( lambda x:  x.primary_key  , 
                             it.ifilter(lambda x:x.__class__==Column,
                                        attrs.values())))[0]):
    attrs['id'] = Column(Integer, primary_key = True)
    print 'no primary key for {0}, making one'.format(name)
  table_type = type(name, (base,), attrs)

  
  return {name:table_type}
