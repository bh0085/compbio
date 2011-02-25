from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import Column, Integer, String, Unicode, DateTime, ForeignKey, UniqueConstraint

def MakeClasses(metadata):
  Base = declarative_base(metadata=metadata)
  class Gbfile(Base):
      __tablename__ = 'gbfile'
      __table_args__ = (UniqueConstraint('name'), {})
      id = Column(Integer, primary_key=True)
      name = Column(String, index = True)
      def __init__(self, name):
          self.name = name
  
  sgjBase = declarative_base(metadata = metadata)
  class Gbfile_GBIDJoin(sgjBase):
    __tablename__ = 'gbfile_gbidjoin'
    __table_args__ = (UniqueConstraint('gbid'), {})

    id = Column(Integer, primary_key=True)
    gbfile = Column(Integer, 
                      ForeignKey('gbfile.id'))
    gbid = Column(Integer,
                  ForeignKey('gbid.id'),
                  index = True)
    def __init__(self, gbfile, gbid):
      self.gbid = gbid
      self.gbfile = gbfile
  return dict(Gbfile = Gbfile, Gbfile_GBIDJoin = Gbfile_GBIDJoin)
  
