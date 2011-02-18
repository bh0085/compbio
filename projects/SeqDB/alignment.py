from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import Column, Integer, String, Unicode, DateTime,ForeignKey,UniqueConstraint

def MakeClasses(metadata):

  Base = declarative_base(metadata=metadata)
  class Alignment(Base):
    __tablename__ = 'alignment'
    __table_args__ = (UniqueConstraint('name'), {})

    id = Column(Integer, primary_key=True)
    source = Column(String)
    name = Column(String)
    def __init__(self, source = 'None'):
        self.source = source

  asjBase = declarative_base(metadata = metadata)
  class Alignment_SequenceJoin(asjBase):
    __tablename__= 'alignment_sequencejoin'
    __table_args__ = (UniqueConstraint('sequence'), {})    
    id = Column(Integer, primary_key=True)
    alignment = Column(Integer, 
                       ForeignKey('alignment.id'))
    sequence = Column(Integer,
                      ForeignKey('sequence.id'))
    def __init__(self, alignment, sequence):
      self.alignment = alignment
      self.sequence = sequence
  return dict(Alignment = Alignment, Alignment_SequenceJoin = Alignment_SequenceJoin)
