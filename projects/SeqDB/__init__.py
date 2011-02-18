from compbio import config
import os 
from sqlalchemy import create_engine
from sqlalchemy.interfaces import PoolListener
from sqlalchemy import MetaData
from sqlalchemy.orm import scoped_session, sessionmaker
import sequence, taxonomy, gbid, alignment
from sqlalchemy.orm.exc import NoResultFound


__all__= ['Session', 'engine','metadata','tryAddUnique']
initiated = False

def setup_db(url = None):
  init(url)
  global metadata
  global engine
  metadata.create_all(bind=engine,checkfirst = False)

def init(url = None):
  if not url:
    url = 'sqlite:///' +  os.path.join(config.root,'sequences/prok/16s.sqlite')
  _initWithURL(url)

def _initWithURL(url):
  global initiated
  if not initiated:
    reload_all()

  class ForeignKeysListener(PoolListener):
    def connect(self, dbapi_con, con_record):
      db_cursor = dbapi_con.execute('pragma foreign_keys=ON')
      
  global engine
  engine = create_engine(url, listeners=[ForeignKeysListener()])
  print 'Binding engine for url ' + url
  global Session
  Session.configure(bind = engine)
  engine = engine
  print Session.connection().engine
  
def MakeMetadata():
  global Session
  global metadata
  global engine
  engine = None
  Session = scoped_session(sessionmaker(autocommit = False))
  metadata =  MetaData()
  global initiated
  initiated = True

#Def reload all classes and put them into the globals
#dict for SeqDB.
#
#Note that this model does not currently support multiple databases at once...
def reload_all():
  MakeMetadata()
  global metadata
  for model in [sequence, gbid, alignment, taxonomy]:
    classes_dict = model.MakeClasses(metadata)
    for k, val in classes_dict.iteritems():
      globals()[k] = val
      if not k in globals()['__all__']:
        globals()['__all__'].append(k)

def tryAddUnique(Obj, **kwargs):
  try: merged = Session.query(Obj.__class__).filter_by(**kwargs).one()
  except NoResultFound:
    merged = Session.merge(Obj)
    Session.flush()
  return merged

