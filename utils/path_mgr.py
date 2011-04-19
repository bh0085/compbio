#!/usr/bin/env python
'''
A few utilities to manage the path stack and to get the filepath of the currently running script

pushd: Push a directory onto the stack.
popd:  Pop a directory from the stack.
getd:  Get the directory of the currently running script.
'''

globals()['dstack'] = []
import inspect, os
def pushd(d = None, f = False):
  '''Push a directory onto the current stack. If no directory is
specified, push the path of the currently runnign script onto
the stack. If f==True, create a directory if none exists.'''
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
  '''Get the directory of the currently running script.'''
  return os.path.dirname(inspect.stack()[1][1])

def popd(f = False):
  '''Pop the topmost directory from the stack.'''
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

