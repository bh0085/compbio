from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import Column, Integer, String, Unicode, DateTime, ForeignKey, UniqueConstraint

def MakeClasses(metadata):
  Base = declarative_base(metadata=metadata)
  class GBID(Base):
      __tablename__ = 'gbid'
      __table_args__ = (UniqueConstraint('name'), {})
      id = Column(Integer, primary_key=True)
      name = Column(String, index = True)
      offset = Column(Integer)
      def __init__(self, name, offset):
          self.name = name
          self.offset = offset
  
  sgjBase = declarative_base(metadata = metadata)
  class GBID_SequenceJoin(sgjBase):
    __tablename__ = 'gbid_sequencejoin'
    __table_args__ = (UniqueConstraint('sequence'), {})

    id = Column(Integer, primary_key=True)
    sequence = Column(Integer, 
                      ForeignKey('sequence.id'),
                      index = True)
    gbid = Column(Integer,
                  ForeignKey('gbid.id'))
    def __init__(self, gbid, sequence):
      self.gbid = gbid
      self.sequence = sequence
  return dict(GBID = GBID, GBID_SequenceJoin = GBID_SequenceJoin)
  
