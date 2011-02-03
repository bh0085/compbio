globals()['dstack'] = []
import inspect
import os
def pushd(d = None, f = False):
  if d == None:
    d = os.path.dirname(inspect.stack()[1][1])
  if d == '':
    d = '.'
  cdir = os.getcwd()
  if not os.path.isdir(d):
    if f:
      os.mkdir(d)
    else:
      raise Exception('Directory {0} does not exist'.format(d))
  os.chdir(d)
  globals()['dstack'].append(cdir)


def getd():
  return os.path.dirname(inspect.stack()[1][1])

def popd(f = False):
  cdir = os.getcwd()
  dname = globals()['dstack'][-1]
  if not os.path.isdir(dname):
    if f:
      os.mkdir(dname)
    else:
      raise Exception('Directory {0} does not exist'.format(dname))
  os.chdir(dname)
  globals()['dstack'].pop()
  return cdir

