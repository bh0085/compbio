#/usr/bin/env python
'''A parser to put the BDNTP virtual embryo files into a python dictionary. Dictionary resembles that from projects.network.io.py >> parseTS().
'''

import os, itertools as it, re
import mpl_toolkits.mplot3d.axes3d as ax3d
import matplotlib.pyplot as plt
from numpy import *
import inspect

def read():
  import inspect
  emb_file = os.path.join(os.path.dirname(inspect.stack()[0][1]), 'embryo.vpc')
  contents = open(emb_file)
  
  rows = []; attrs = {}
  for l in contents.xreadlines():
    if l[0:2] == '##': continue
    if l[0] == '#':
      k, v = [x.strip() for x in l[2:].split('=')]
      attrs[k] = v
      continue
    rows.append(l)

    
  contents.close()
  cinfo = attrs['column_info']
  cinfo_names = ['cname', 'type', 'short_name', 'long_name', 'data_type','subc']
  cinfo = [[x.replace('"','') for x in row.split(',')]
           for row in  cinfo[1:-1].split(';')]
  cinfo_out = dict([(celt[0].strip(), dict(zip(cinfo_names, (celt)))) for celt in cinfo])
  attrs['column_info'] = cinfo_out
  attrs['column'] = attrs['column'][1:-1].split(',')

  col_inames = attrs['column_info'].keys()

  attrs['column'] = [c.replace('"', '').strip() for c in attrs['column']]
  g_colnames =  attrs['column']

  scols =  sorted(enumerate(g_colnames), key = lambda x: x[1])
  col_grps = dict([((k,list(g)))
              for k,g in it.groupby(
      scols,
      lambda x:re.compile('[^_]*').search(x[1]).group())
              ])
  
  
  
  column_gnames = set(col_inames)
  column_miscnames = set( col_grps.keys()).difference(column_gnames)

  gene_cols =dict( [(k, dict(
        info = attrs['column_info'][k],
        idxs = array(zip(*col_grps[k])[0], int),
        fullnames = zip(*col_grps[k])[1],
        steps = [n[-1] for n in zip(*col_grps[k])[1]]))
               for k in column_gnames])


  misc_cols = dict([(k,dict(idxs =array( zip(*col_grps[k])[0], int),
                            fullnames = zip(*col_grps[k])[1],
                            steps = [n[-1] for n in zip(*col_grps[k])[1]]))
                    for k in column_miscnames])
  

  max_std_col = max(it.chain(*[x['idxs'] 
                               for x in list(it.chain(misc_cols.values(),
                                                      gene_cols.values()))]))
  row_nns = [x.split(',')[max_std_col + 1:] for x in rows]
  rows = array([array(x.split(',')[0:max_std_col+1], float) for x in rows])

  return gene_cols, misc_cols, rows, row_nns

def show(rows):

  fig = plt.figure(1)
  try:plt.clf()
  except Exception, e: print e
  
 
  xs = rows[:,misc_cols['x']['idxs'][0]]
  ys = rows[:,misc_cols['y']['idxs'][0]]
  zs = rows[:,misc_cols['z']['idxs'][0]]
  ax = fig.add_subplot('111', projection='3d')

  #MAKE A SINGLE TIMEPOINT OR A TIME SERIES
  make_gif = True
  if make_gif: rng = range(6)
  else: rng = [0]
  nsteps = 50
  stepwid = len(rows) / nsteps * 2

  for step in [-2]: #range(nsteps):
    start, stop = step * len(rows) / (nsteps-1) -stepwid, step * len(rows) / (nsteps-1) +stepwid
    start, stop = max([0,start]), min([stop,len(rows)])
    xinds = argsort(xs)[start: stop]

    tsums = rows[:,gene_cols['twi']['idxs'][0]]*50.
    ssums = rows[:,gene_cols['sna']['idxs'][0]]*50.
    ax.scatter(xs[xinds], ys[xinds], zs[xinds],
               s = ssums[xinds])
    ax.scatter(xs[xinds] + 4, ys[xinds], zs[xinds],
               s = tsums[xinds] , color = 'red' )

    #ax.set_xlim(min(xs), max(xs))
    ax.set_ylim3d(min(ys),max(ys))
    ax.set_xlim3d(min(xs),max(xs))
    ax.set_zlim3d(min(zs),max(zs))
    plt.savefig(os.path.join(os.path.dirname(inspect.stack()[0][1]),'twi_sna_{0}'.format('%03i' % (step,))))
    plt.savefig(os.path.join(os.path.dirname(inspect.stack()[0][1]),'twi_sna_{0}'.format('%03i' %( 2*nsteps - step))))
    ax.clear()
  return rows, attrs


def make_movie():
  #
  # Now that we have graphed images of the dataset, we will stitch them
  # together using Mencoder to create a movie.  Each image will become
  # a single frame in the movie.
  #
  # We want to use Python to make what would normally be a command line
  # call to Mencoder.  Specifically, the command line call we want to
  # emulate is (without the initial '#'):
  # mencoder mf://*.png -mf type=png:w=400:h=300:fps=25 -ovc lavc -lavcopts vcodec=mpeg4 -oac copy -o out.avi
  # See the MPlayer and Mencoder documentation for details.
  #

  filepath = os.path.dirname(inspect.stack()[0][1])
  command = ('mencoder',
             'mf://{0}/*.png'.format(filepath),
             '-mf',
             'type=png:w=800:h=600:fps=30',
             '-ovc',
             'lavc',
             '-lavcopts',
             'vcodec=mpeg4',
             '-oac',
             'copy',
             '-o',
             'out.avi')

  import subprocess
  #os.spawnvp(os.P_WAIT, 'mencoder', command)

  print "\n\nabout to execute:\n%s\n\n" % ' '.join(command)
  subprocess.check_call(command)

  print "\n\n The movie was written to 'output.avi'"
  
  print "\n\n You may want to delete *.png now.\n\n"
