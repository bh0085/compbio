import os
import socket
root = os.environ['COMPBIO_PATH']

def absPath(localPath):
  global root
  return os.path.join(root, localPath)

def dataPath(url):
  if url.count(':') == 0:
    host_name = 'localhost'
    volume_name = 'cb'
    localpath = url
  elif url.count(':') == 2:
    host_name = url.split(':')[0]
    volume_name = url.split(':')[1]
    localpath = url.split(':')[2]
  
  if volume_name == '/':volume_prefix = '/'
  elif volume_name == 'cb':volume_prefix = globals()['root']
  elif volume_name == '~':volum_prefix = os.environ['HOME']
  else: volume_prefix = os.path.join('/Volumes', volume_name)
  assert host_name == socket.gethostname() or \
      host_name == 'localhost'
  return os.path.join(os.path.join(volume_prefix,'data'), localpath)

def dataURL(localpath, volume_name = 'cb',  host_name = 'localhost'):
  return ':'.join([host_name, volume_name, localpath]) 

def sqlitePath(dbname, **kwargs):
  url = dataURL(os.path.join('dbs',dbname+'.sqlite'),**kwargs)
  return 'sqlite:///'+dataPath(url)
def postgresDefault(host = None):
  if host == 'broad':
    return 'postgres://benh@node1386.broadinstitute.org/benh'
  else:
    return 'postgres://bh0085@localhost/bh0085'
def postgresPath(dbname, host = None):
  if host == 'broad':
    return 'postgres://benh@node1386.broadinstitute.org/'+dbname
  else:
    return 'postgres://bh0085@localhost/'+dbname