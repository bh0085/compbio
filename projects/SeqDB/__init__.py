from compbio import config
import os 
from sqlalchemy import create_engine
from sqlalchemy.interfaces import PoolListener
from sqlalchemy import MetaData
from sqlalchemy.orm import scoped_session, sessionmaker
import sequence, taxonomy, gbid, alignment, gbfile
from sqlalchemy.orm.exc import NoResultFound

def dbname(name):
  name_str = 'db_' + name
  if not name_str in globals().keys():
    raise Exception( 'db not yet initialized for name ' + name)
  return globals()[name_str]


def _register_dbconnected(db, name):
  globals()['db_' + name] = db


class SeqDBI():
  __all__= ['Session', 'engine','metadata','tryAddUnique']

  def setup_db(self, url = None):
    self.metadata.create_all(bind=self.engine,checkfirst = False)

  def __init__(self, url = None):
    if not url:
      url = 'sqlite:///' +  os.path.join(config.root,'sequences/prok/16s.sqlite')
    self.name = os.path.splitext(os.path.basename(url))[0]
    self.reload_all()
    _register_dbconnected(self, self.name)

    class ForeignKeysListener(PoolListener):
      def connect(self, dbapi_con, con_record):
        db_cursor = dbapi_con.execute('pragma foreign_keys=ON')      
    self.engine = create_engine(url, listeners=[ForeignKeysListener()])
    print 'Binding engine for url ' + url
    self.Session = scoped_session(sessionmaker(autocommit = False, 
                                               bind = self.engine))
    self.url = url
  
  def reload_all(self):
    self.metadata = MetaData()
    for model in [sequence, gbid, alignment, taxonomy, gbfile]:
      classes_dict = model.MakeClasses(self.metadata)
      for k, val in classes_dict.iteritems():
        self.__dict__[k] = val

  def tryAddUnique(self,Obj, **kwargs):
    try: merged = self.Session.query(Obj.__class__).filter_by(**kwargs).one()
    except NoResultFound:
      merged = self.Session.merge(Obj)
      self.Session.flush()
    return merged

