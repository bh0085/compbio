import os, pipes, socket, subprocess as spc
root = os.environ['COMPBIO_PATH']
import compbio.utils.remote_utils 

def remotePath(abspath, host = 'tin', root = 'comp'):
  '''Get the location of the a file on the remote file system
from one of the roots ['compbio', 'programming']'''
  
  if root == 'prog':
    subfun = progPath
    subterm =  'progPath'
  else:
    subfun = compPath
    subterm =  'compPath'    
  
  if host == None:
    #IS THIS WHY THE LOCAL CALLS TO BSUB ARE FAILING?
    print abspath
    print subfun(abspath)
    return subfun(abspath).strip()

  scr = pipes.quote('''
echo `python -c {0}`'''.format(pipes.quote('''
import compbio.config as config
import os, inspect
print config.{1}('{0}', absolute = True)
'''.format( subfun(abspath),
            subterm
            ))))

  
  ssh_scr = 'ssh {1} {0}'.format(scr, host)
  out = spc.Popen(ssh_scr, shell = True, stdout = spc.PIPE).\
      communicate()[0]
  return out.strip()

def absPath(localPath):
  global root
  return os.path.join(root, localPath)

def compPath(path, absolute = False):
  if absolute:
    return os.path.join(os.environ['COMPBIO_PATH'], path)
  else:
    return os.path.relpath(os.path.abspath(path), \
                             os.environ['COMPBIO_PATH']) 
  
def progPath(path, absolute = False):
  if absolute:
    return os.path.join(os.environ['PROGRAMMING_PATH'], path)
  else:
    return os.path.relpath(os.path.abspath(path), \
                             os.environ['PROGRAMMING_PATH']) 
  
def getTempPath():
  return dataPath('temp')

def dataPath(url, make = True):
  if url.count(':') == 0:
    host_name = 'localhost'
    volume_name = 'cb'
    localpath = url
  elif url.count(':') == 2:
    host_name = url.split(':')[0]
    volume_name = url.split(':')[1]
    localpath = url.split(':')[2]

  if host_name == '':
    host_name = 'localhost'
  if volume_name == '':
    volume_name = 'cb'

  if volume_name == '/':volume_prefix = '/'
  elif volume_name == 'cb':volume_prefix = globals()['root']
  elif volume_name == '~':volum_prefix = os.environ['HOME']
  else: volume_prefix = os.path.join('/Volumes', volume_name)
  if host_name == socket.gethostname() or \
        host_name == 'localhost':
    path = os.path.join(os.path.join(volume_prefix,'data'), localpath)
    if make and not os.path.isdir(os.path.dirname(path)):
      os.makedirs(os.path.dirname(path))
  else:
    rutils =  compbio.utils.remote_utils 
    path = rutils.remote_datapath(localpath, 
                                  host_name, 
                                  volume = volume_name)
  return path

def dataURL(localpath, volume_name = 'cb',  host = 'localhost'):
  return ':'.join([host, volume_name, localpath]) 

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

def scriptInputPath(name):
  return os.path.join(os.path.join(globals()['root'], 
                                   'scripts/scr_inputs'),
                      name)
def scriptOutputPath(name):
  return os.path.join(os.path.join(globals()['root'],
                                   'scripts/scr_outputs'),
                      name)
