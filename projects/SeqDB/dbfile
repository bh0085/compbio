from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import Column, Integer, String, Unicode, DateTime, ForeignKey, UniqueConstraint

def MakeClasses(metadata):
  Base = declarative_base(metadata=metadata)
  class File(Base):
      __tablename__ = 'file'
      __table_args__ = (UniqueConstraint('name'), {})
      id = Column(Integer, primary_key=True)
      name = Column(String, index = True)
      def __init__(self, name):
          self.name = name
  
  sgjBase = declarative_base(metadata = metadata)
  class File_GBIDJoin(sgjBase):
    __tablename__ = 'gbid_sequencejoin'

    id = Column(Integer, primary_key=True)
    file = Column(Integer, 
                      ForeignKey('file.id'),
                      index = True)
    gbid = Column(Integer,
                  ForeignKey('gbid.id'))
    def __init__(self, file, gbid):
      self.gbid = gbid
      self.file = file
  return dict(File = File, File_GBIDJoin = File_GBIDJoin)
  
