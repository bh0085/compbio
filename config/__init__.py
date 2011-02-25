import os
import socket
root = os.environ['COMPBIO_PATH']

def absPath(localPath):
  global root
  return os.path.join(root, localPath)

def dataPath(url):
  if url.count(':') == 0:
    host_name = socket.gethostname()
    volume_name = '/'
    localpath = url
  elif url.count(':') == 2:
    host_name = url.split(':')[0]
    volume_name = url.split(':')[1]
    localpath = url.split(':')[2]
  
  if volume_name != '' and volume_name != '/':
    volume_prefix = os.path.join('/Volumes/',volume_name)
  else:
    volume_prefix = '/'

  assert host_name == socket.gethostname() or \
      host_name == 'localhost'

  return os.path.join(os.path.join(volume_prefix,'data'), localpath)

def dataURL(localpath, volume_name = '/',  host_name = None):
  if host_name == None:
    host_name = socket.gethostname()
  assert host_name == socket.gethostname()
  return ':'.join([host_name, volume_name, localpath]) 

def sqlitePath(dbname):
  return 'sqlite:////'+os.path.join(globals()['root'], 'data/dbs/'+dbname+'.sqlite')
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
