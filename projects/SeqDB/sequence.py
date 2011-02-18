from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import Column, Integer, String, Unicode, DateTime, ForeignKey, UniqueConstraint

def MakeClasses(metadata):
  Base = declarative_base(metadata=metadata)
  class Sequence(Base):
    __tablename__ = 'sequence'
    __table_args__ = (UniqueConstraint('name'), {})

    id = Column(Integer, primary_key=True)
    sequence = Column(String, nullable = False)
    name = Column(String, nullable = False)
    source = Column(String)
    description = Column(String)
    def __init__(self, sequence,name,source, description):
      self.sequence = sequence
      self.name = name
      self.source = source
      self.description = description
  return dict(Sequence = Sequence)
