import matplotlib.pyplot as plt

def run_subroutine():
  print 'this subroutine is running'
  return 0

def run():
  print 'this program is running'
  out = run_subroutine()
  print 'subroutine returned {0}'.format(out)

  raise Exception()

  print 'hi'

  return

run()
