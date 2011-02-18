__all__ = ['Session', 'engine', 'metadata']
engine = None
Session = scoped_session(sessionmaker(autocommit = True))
metadata =MetaData()
