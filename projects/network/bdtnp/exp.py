#!/usr/bin/env python
'''
A few random experiments for the BDTNP data to get a feel for
what's going on.

p_m_correlation:   plots the correlations between protein and mRNA
  levels for the four genes for which it is available and saves
  a figure for later reference.

show_multi 
  Shows the expression levels for a bunch of different mRNAs in
  the form of a 3 dimensional scatterplot.

show_3d
  Shows a single gene in a 3 dimensional scatterplot over the 
  embryo.

cluster_tissues
  Representing nuclei by the expression levels of a set d of genes
  at a particular time point, clusters nuclei into medioid clusters
  using the max-sum algorithm from gifford's class.

  Uses the online interface for now.

'''

from compbio.projects.network import io as nio, utils as nu
from compbio.utils import colors as mycolors, path_mgr as pm
from compbio.utils import memo as mem, seismic as seismic
import compbio.config as cfg
import matplotlib.pyplot as plt
import pickle
from numpy import *
import numpy as np, textwrap as tw
import cs874.bsub_clusters as bcl
import os, inspect, itertools as it
import scipy.io as sio
import scipy.stats as ss

'''
Short documentation:

ll = launchClusters()
ll.launch()
ll.fetch_start()
all_members,ct_data = cluster_exprs()
cluster_exprs(all_members, ct_data)

'''

if not os.path.isdir(cfg.dataPath('figs/bdtnp')):
  os.mkdir(cfg.dataPath('figs/bdtnp'))


def getClusteringInputs( ncluster =1000, host = 'tin', 
                           reset = 0, step = 10, exemp_time = 'all',
                           doplot = False):
  mrnas = nio.getBDTNP()
  misc = nio.getBDTNP(misc = True)
  
  vals = array([v['vals'] for v in mrnas.values()])
  gvars = var(vals, 1)
  gminvars = np.min(gvars,1)
  gmedvars = median(gvars,1)

  min20 = argsort(gminvars)[::-1][:20]
  med20 = argsort(gmedvars)[::-1][:20]

  int20 = set(min20).intersection(set(med20))
  xgenes = array(list(int20))

  cell_data = vals[xgenes].transpose(1,2,0)
  scd = shape(cell_data)
  #times = reshape(zeros(shape(cell_data[0:2]))[:,:,newaxis , arange(shape(cell_data[1]))
  #                    , (prod(shape(cell_data)[0:2])))
  xycoords = (arange(scd[0])[:,newaxis,newaxis]*[1,0] +\
                arange(scd[1])[newaxis,:,newaxis]*[0,1])
  cell_data = reshape(cell_data, (prod(shape(cell_data)[0:2]), shape(cell_data)[2] ))
  xy_data = reshape(xycoords, (prod(scd[0:2]),2 ))
    
  if exemp_time == 'all':
    inds = arange(len(cell_data))
  else:
    inds = arange(len(cell_data))[nonzero(equal(xy_data[:,1],exemp_time))[0]]
  
  np.random.seed(1)
  np.random.shuffle(inds)
  rand_thousand = inds[0:ncluster]
  
  sim_data = cell_data[rand_thousand]
  sim_xy = xy_data[rand_thousand]
  t = [ mean(sim_data, 0), std(sim_data,0)]
  t[1][equal(t[1],0)] = 0
  metric = 'neg_dist'
  sims = similarity(sim_data, transform = t, method = metric)

  name = 'll_{0}_{1}_{2}'.format(metric,ncluster,exemp_time)

  d_in = []
  percs = logspace(.1,1.5,8)
  for p in percs:
    d_in.append(dict(similarities = sims,
                     self_similarity = ss.scoreatpercentile(sims, p),
                     metric = metric
                     ))
  return d_in


def launchClusters(ncluster = 2000):
  dicts = list(it.chain(*[getClusteringInputs(ncluster = ncluster, host = 'tin',exemp_time = t)
                          for t in [0,1,2,3,4,5,'all']]))
  launcher = bcl.launcher(dicts, host = 'tin', name = 'apclustering_n{0}'.format(ncluster))
  launcher.launch()
  return launcher


def p_m_correlation():
  prots = nio.getBDTNP(protein = True)
  mrnas = nio.getBDTNP()
  
  matched = set(mrnas.keys()).intersection(set(prots.keys()))
  pairs = [(prots[k] , mrnas[k], k) for k in matched]



  f = plt.figure(0)
  f.clear()
  f.suptitle('mRNA and Protein Levels from BDTNP at six times in ~6000 cells', fontsize = 22)
  nx = ny = ceil(sqrt(len(pairs)))
  
  shp = shape(mrnas.values()[0]['vals'])
  colors = mycolors.getct(shp[1])
  shr = None
  for i, p in enumerate(pairs):
    
    ax = f.add_subplot('{0:g}{1:g}{2:g}'.format(nx, ny , i+1),
                       sharex = shr,sharey = shr)
    if not shr: shr = ax
    fbid = p[-1]
    #ax.set_title('{2}'.format(\
    #    fbid, nu.gene_symbol(fbid), tw.fill(nu.gene_biology(fbid), 75)))
    ax.grid(True, alpha = .2)
    ax.annotate(nu.gene_symbol(fbid),xy = [.02,.98], 
                xycoords = 'axes fraction', size = 25, va = 'top')
    mu = corrcoef(p[0]['vals'][::,:].flatten(),p[1]['vals'][::,:].flatten())
    ax.annotate('$\mu = {0:.2g}$'.format(mu[0,1]),xy = [.98,.98],
                xycoords = 'axes fraction', size = 25,ha = 'right', va = 'top')

    if mod(i, nx) >0: 
      plt.setp( ax.get_yticklabels(), visible=False)
    else:  ax.set_ylabel('mrna expression level')
      #plt.setp( ax.get_ylabel(), visible=False)
    if floor(i/nx) < (ny -1) : 
      plt.setp( ax.get_xticklabels(), visible=False)
    else:  ax.set_xlabel('protein expression level')

      #plt.setp( ax.get_xlabel(), visible=False)

    for j in range(shp[1]):
      ax.scatter(p[0]['vals'][::,j],p[1]['vals'][::,j],
                 s = 20,alpha = .2,color = colors[j])
      
  f.savefig(cfg.dataPath('figs/network/mrna_protein_levels.tiff',
                            ),format = 'tiff')
    
                       


def recall_c2(**kwargs):
  '''
A kludgy wrapper to store the clustering results for later
without modifying the original mess of a program, c2...
'''
  def setC2(**kwargs):
    ll = c2(**mem.sr(kwargs))
    result =  c2(ll, **mem.sr(kwargs))
    return result
  return mem.getOrSet(setC2, 
                      **mem.rc(kwargs,
                               name = 'default_c2_settings',
                               on_fail = 'compute'))
def c2( launcher = None, ncluster =2000, host = 'tin', 
        reset = 0, step = 10, exemp_time = 'all',
        doplot = False ,**kwargs):
  mrnas = nio.getBDTNP()
  misc = nio.getBDTNP(misc = True)
  
  vals = array([v['vals'] for v in mrnas.values()])
  gvars = var(vals, 1)
  gminvars = np.min(gvars,1)
  gmedvars = median(gvars,1)

  min20 = argsort(gminvars)[::-1][:20]
  med20 = argsort(gmedvars)[::-1][:20]

  int20 = set(min20).intersection(set(med20))
  xgenes = array(list(int20))

  cell_data = vals[xgenes].transpose(1,2,0)
  scd = shape(cell_data)
  #times = reshape(zeros(shape(cell_data[0:2]))[:,:,newaxis , arange(shape(cell_data[1]))
  #                    , (prod(shape(cell_data)[0:2])))
  xycoords = (arange(scd[0])[:,newaxis,newaxis]*[1,0] +\
                arange(scd[1])[newaxis,:,newaxis]*[0,1])
  cell_data = reshape(cell_data, (prod(shape(cell_data)[0:2]), shape(cell_data)[2] ))
  xy_data = reshape(xycoords, (prod(scd[0:2]),2 ))
    
  if exemp_time == 'all':
    inds = arange(len(cell_data))
  else:
    inds = arange(len(cell_data))[nonzero(equal(xy_data[:,1],exemp_time))[0]]
  
  np.random.seed(1)
  np.random.shuffle(inds)
  rand_thousand = inds[0:ncluster]
  
  sim_data = cell_data[rand_thousand]
  sim_xy = xy_data[rand_thousand]
  t = [ mean(sim_data, 0), std(sim_data,0)]
  t[1][equal(t[1],0)] = 0
  metric = 'neg_dist'
  sims = similarity(sim_data, transform = t, method = metric)

  name = 'll_{0}_{1}_{2}'.format(metric,ncluster,exemp_time)
  def setLauncher(**kwargs):
    sims= kwargs.get('sims')
    metric = kwargs.get('metric')
    name = kwargs.get('name')
    d_in = []
    percs = logspace(.1,1.5,8)
    for p in percs:
      d_in.append(dict(similarities = sims,
                       self_similarity = ss.scoreatpercentile(sims, p),
                       metric = metric
                       ))

    launcher = bcl.launcher(d_in, host = host, name = name)
    return launcher  
  if launcher == None:
    output = mem.getOrSet(setLauncher,
                          **mem.rc(dict(sims = sims, metric = metric,
                                        name = name,
                                        hardcopy = True,
                                        reset = reset,
                                        hard_reset = False,)))  
    return output



  def setC2(launcher = launcher, **kwargs):
    if launcher == None:
      raise Exception()
    else:
      output = launcher.output()
    return output
    #It appears that the bsub process failed for the first output.
    #No big deal. Debug later.
  
  output = mem.getOrSet(setC2,
                        **mem.rc(dict(harcopy = True,
                                      launcher = launcher,
                                      reset = reset,
                                      on_fail = 'compute',
                                      hard_reset = False,
                                      name =  'c2'+ name )))
  all_inds = array([  squeeze(o['inds']) for o in output[:] ])
  

  xs = misc['x']['vals'][zip(*xy_data)] #zip(*sim_xy)]
  ys = misc['y']['vals'][zip(*xy_data)] #zip(*sim_xy)]
  zs = misc['z']['vals'][zip(*xy_data)] #zip(*sim_xy)]
  
  colors =array( mycolors.getct(shape(all_inds)[1]) )
  f = plt.figure(0)
  f.clear()
  
  all_tps = range(scd[1])
  nc = len(all_inds)
  nt = len(all_tps)

  all_members = []
  for i, inds in enumerate(all_inds):
    #compute similarity matrices 1000 at a time:
    exemplars = sim_data[list(set(list(inds)))]
    sim = similarity(cell_data, 
                   exemplars, 
                   transform = t,
                   method = metric)
    closest = argmax(sim, 1)
    all_members.append(closest)
    
    
    if doplot:
      for j, tp in enumerate(all_tps):
        ax = f.add_axes( [float(j)/nt,float(i) /nc,1./nt, 1. /nc] )
        ax.set_yticks([])
        ax.set_xticks([])
        i_sub = nonzero(equal(xy_data[:,1], j) * greater(ys,0))[0]
        cs = colors[closest[i_sub]]
        x = xs[i_sub]
        z = zs[i_sub]
        plt.scatter(x[::step],z[::step], 40,alpha = .75, c = cs[::step], edgecolor = 'none')
    
  ct_data = xy_data
  return all_members, ct_data

def cluster_exprs(all_members, ct_data,
                  do_plot = False,
                  cluster_type = '4d',
                  cluster_id = 4):
  mrnas = nio.getBDTNP()
  misc = nio.getBDTNP(misc = True)

  c = all_members[cluster_id]
  c_unq = set(list(c))
  

  tissues = dict([('t_{0}'.format(i) , dict(cts = ct_data[equal(c,elt)]))
                  for i, elt in enumerate(c_unq)])
  
  nt = 6
  counts = array([[sum(equal(v['cts'][:,1],t))
                   for t in range(nt) ] 
                  for v in tissues.values() ])
  

  if do_plot:
    f = plt.figure(1)
    f.clear()
  
    ax1 = f.add_subplot('121')
    ax2 = f.add_subplot('122')
    seismic.seismic(counts , ax = ax1,stacked = True,colors = mycolors.getct(len(counts)))
    #seismic.seismic(np.sort(counts,0) , ax = ax2,stacked = False,colors = mycolors.getct(len(counts)))
    ax2.hist(np.sum(counts,1))
    
  
  all_exprs = {}
  for t, v in tissues.iteritems():
    ct_all = v['cts']
    
    for time in set([c[1] for c in ct_all]):
      ct = [ct for ct in ct_all if ct[1] == time]

      exprs =dict( [(k,elt['vals'][zip(*ct)]) for k, elt in mrnas.iteritems()])
      ys = misc['y']['vals'][zip(*ct)] #zip(*sim_xy)]
      zs = misc['z']['vals'][zip(*ct)] #zip(*sim_xy)]
      xs = misc['x']['vals'][zip(*ct)] #zip(*sim_xy)]

    
      f = plt.figure(1)
      f.clear()
      ax1 = f.add_subplot('121', title = 'X-Z axis view for tissue {0}'.\
                            format(t))
      ax2 = f.add_subplot('122',title = 'Y-Z axis view for tissue {0}'.\
                            format(t))
      ax1.scatter(xs, zs)
      ax2.scatter(ys, zs)
      
      v['exprs'] = exprs
      all_exprs['tiss_{0}_time_{1}'.format(t,time)]=exprs
      
      sio.savemat(open(cfg.dataPath('soheil/expression_c{0}_n{1}_tissue{2}_time{3}.mat'.\
                                      format(cluster_type,cluster_id,t,time)),'w'),
                  exprs)
      f.savefig(open(cfg.dataPath('soheil/expression_c{0}_n{1}_tissue{2}_time{3}.tiff'.\
                                      format(cluster_type,cluster_id,t,time)),'w'))
    
    
      
  exprs_out = dict([( k, [ mean(sub[k]) for sub in all_exprs[k].values() ]) 
                    for k in all_exprs.keys() ])

  sio.savemat(open(cfg.dataPath('soheil/expression_c{0}_n{1}_intercluster.mat'.\
                                    format(cluster_type,cluster_id)),'w'),
              exprs_out)
  
  raise Exception()


def cluster_tissues(nx = 20,ny = 500, timepoint = -1,
                    step = 4,
                    sim = 'neg_dist', 
                    imshow_sims = False,
                    scatter_sims = False,
                    hist_sims = False,
                    do_cluster= True,
                    do_show = True, cstep = -1):
  '''Cluster ny nuclei by the values of the nx mRNAs with highest
  variance. Uses the medioids method with number of clusters set
  by exemplar self simalarity as outlined in 6.874 and implemented
  at http://www.psi.toronto.edu/affinitypropagation 


  imaging:
  imshow_sims
  scatter_sims
  hist_sims
  do_show

  numerics:
  nx:           number of genes to cluster upon
  ny:           number of cells in the clusterin
  timepoint:    which time to use for cluster computation
  step:         how many genes to skip when showing results
  
  
  So far I have implemented a distance based similarity and a 
  '''
  mrnas = nio.getBDTNP()
  misc = nio.getBDTNP(misc = True)
  shp = shape(mrnas.values()[0]['vals'])

  #choose to look only at one timepoint
  stds = [std(m['vals'][:,timepoint]) for m in mrnas.values()]
  vsort = argsort(stds)[::-1]
  xinds = vsort[:nx]

  #Choose the most variable factors and use them as the 
  #underlying variables from which to construct a similarity
  nuclei =array([ mrnas.values()[idx]['vals'][:,timepoint]
                  for idx in xinds]).T

  t = [ mean(nuclei, 0), std(nuclei,0)]
  t[1][equal(t[1],0)] = 0
  sims = similarity(nuclei, transform = t, method = sim)
  cluster_inds = array(floor(linspace(0,len(nuclei)-1, ny)), int)  
  cluster_training = sims[cluster_inds,:][:,cluster_inds]

  f = plt.figure(0)

  #, projection = '3d')

  if scatter_sims:
    ax = f.add_subplot(111)
    scatterx = [cluster_sims[i] for i in range(ny) for j in range(ny)]
    scattery = [cluster_sims[j] for i in range(ny) for j in range(ny)]
    ax.scatter(scatterx, scattery, s =3, alpha = .1)

  if imshow_sims:
    ax = f.add_subplot(111)
    cmap = mycolors.blackbody()
    ax.imshow(cluster_sims, cmap = cmap, interpolation = 'nearest')

  if hist_sims:
    ax = f.add_subplot(111)
    csf = cluster_sims.flatten()
    csf -= max(csf)
    csf *= -1
    h = histogram(log10(1+csf), bins = 100)
    ax.plot(h[1][:-1],h[0])


  cluster(cluster_training, ss.scoreatpercentile(cluster_training,.2) )
  
  fopen = open(cfg.dataPath('bdtnp/clustering/nuclei/idxs'))
  lines = fopen.readlines()
  c = [int(l.strip()) for l in lines]
  c_training_exemplars = set(c)
  exemplar_inds = [cluster_inds[i] for i in c_training_exemplars]
  #I am being a bit lazy with subscripting here because I just assume
  #that the similarity is symmetric... I suppose I could let it be 
  #asymmetric if I liked

  
  exemplars = nuclei[exemplar_inds,:]
  all_sims = similarity(nuclei,  exemplars,
                        transform = t, 
                        transform_exemplars = True,
                        method = sim)
  assignments = np.argmax(all_sims,1)


  ne = len(c_training_exemplars)
  colors = array(mycolors.getct(len(c)))
  colors = array(colors)


  if do_show:
    for tp in range(shape(mrnas.values()[0]['vals'])[1])[-1:]:
      try: f.clear()
      except Exception, e: print 'Weird 3d plotting error. Alas'
      nuclei =array([ mrnas.values()[idx]['vals'][:,tp]
                      for idx in xinds]).T
      all_sims = similarity(nuclei,  exemplars,
                            transform = t, transform_exemplars = True,
                            method = sim)
      assignments = np.argmax(all_sims,1)


      ax = f.add_subplot(111)
      #colors = [colors[i] for i in c]
      xs = misc['x']['vals'][::step,0]
      ys = misc['y']['vals'][::step,0]
      zs = misc['z']['vals'][::step,0]
      ax.scatter(xs, zs,s= 50, color =colors[assignments[::step]])
      #ax.set_title('''Virtual embryo cell (n={2}) clusters 
#from similarities derived from {0} genes. 
#Clusters derived at T = {1}, shown at T = {3}.'''\
#                     .format(nx,timepoint, len(xs),tp))
    
      f.savefig(cfg.dataPath('figs/bdtnp/cluster_movie{0:02d}.tiff'.format(tp)), format = 'tiff')
      
      
def similarity(query,exemplars = None, 
               transform = (0,1),  transform_exemplars = True, 
               method = 'neg_dist'):
  '''Compute the similarities of matrix query to matrix exemplars. If
exemplars is unspecified, compute the similarity of matrix query to 
itself. Query should be an [m, nx] matrix and exemplars should be an [n, nx]
matrix.

args:
  query:     [n, nx]  matrix for which to compute similarities
  exemplars: [m, nx]  exemplar matrix 
  transform:  [t,scl]  transform to apply to xcoords of query to match exemplars
                      t is subtracted from query
                      scl is then divided by
  transform_exemplars [True/False] transform exemplars or only query.

  method:    ['neg_diff','angle']
  

'''
  
  qscl = (query - transform[0]) / transform[1]
  if exemplars == None: 
    escl = array(qscl)
  else:
    if transform_exemplars:  escl = (array(exemplars) - transform[0]) / transform[1]
    else: escl = array(exemplars)
  

  #make sure that we are not allowing all zero columns into either 
  #our query or exemplar vectors in the similarity computation
  good = nonzero(greater(var(qscl, 0) * var(escl,0), 0))[0]

  print 'Print used {0} genes out of {1}'.format(len(good), shape(query)[1])
  #compute a similarity
  if method == 'angle': 
    cq_units = qscl[:,good] / sqrt(sum(qscl[:,good]**2, 1))[:,newaxis]
    ce_units = exemplars[:,good] / sqrt(sum(exemplars[:,good]**2, 1))[:,newaxis]
    sims= dot(cq_units, ce_units.T)
  elif method == 'neg_dist':
    cq_norms = sum(qscl[:,good]**2,1)
    ce_norms = sum(escl[:,good]**2,1)
    cluster_dists = cq_norms[:,newaxis] + ce_norms[newaxis,:] \
        - 2 * dot(qscl[:,good], escl[:,good].T)
    sims = - cluster_dists
  
  return sims

def make_movie(moviedir = cfg.dataPath('figs/bdtnp')):
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

  import subprocess

  filepath = os.path.dirname(inspect.stack()[0][1])
  command = ('mencoder',
             'mf://{0}/*.png'.format(moviedir),
             '-mf',
             'type=png:w=800:h=600:fps=2',
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
  


def cluster(similarities, self_sim):
  if not os.path.isdir( cfg.dataPath('bdtnp/clustering/nuclei/')):
    os.mkdir( cfg.dataPath('bdtnp/clustering/nuclei/'))

  ny = len(similarities)
  simfile = open(\
    cfg.dataPath('bdtnp/clustering/nuclei/Similarities.txt'),'w')
  ssfile = open(\
    cfg.dataPath('bdtnp/clustering/nuclei/Preferences.txt'),'w')
  simlines = ['{0:05d}   {1:05d}  {2:g}\n'.\
                format(i+1, j+1, similarities[i,j]) 
              for i in range(ny) for j in range(ny) if i != j]
  for s in simlines: simfile.write(s)
  preflines = ['{0:0.8g}\n'.format(self_sim) for i in range(ny)]
  for p in preflines: ssfile.write(p)
  
  ssfile.close()
  simfile.close()

def show_multi(timepoint = -1):
  mrnas = nio.getBDTNP()
  misc = nio.getBDTNP(misc = True)
  shp = shape(mrnas.values()[0]['vals'])

  #choose to look only at one timepoint
  stds = [std(m['vals'][:,timepoint]) for m in mrnas.values()]
  
  f = plt.figure(0)
  try: f.clear()
  except Exception, e: print 'hi'
  ax = f.add_subplot(111, projection = '3d')
  vsort = argsort(stds)[::-1]
  
  n = 10
  colors = mycolors.getct(n)
  for i in arange(n):
    step = argmax(np.sum(mrnas.values()[vsort[i]]['vals'],0))
    show_3d(mrnas.keys()[vsort[i]], 
            step = step, skip = 20, ax = ax, ofs =10*random.rand(3),
            color = colors[i])

def show_3d(fbid, scatter = False,
            step = 0, skip = 10,ofs = [0,0,0], ax = None, **kwargs):
  '''Plot the expression of a given mrna in 3d.
If you pass it an axes, you can get hte same axis back and run the script again to plot other genes.

step, skip set which timestep and how many nuclei to plot.

ax sets the plotting ax
**kwargs sets keywords in the scatterplot
'''

  mrnas = nio.getBDTNP()
  misc = nio.getBDTNP(misc = True)  
  
  assert fbid in mrnas.keys(), 'No BDTNP data for {0}'.format(fbid)

  xs,ys,zs = array([misc[k]['vals'][::skip,step] 
                    for k in ('x','y','z')]) + ofs[:,newaxis]

  sizes = array(mrnas[fbid]['vals'][::skip,step])
  sizes -= min(sizes)
  sizes/= max(sizes)/10. if max(sizes)  > 0 else 1
  
  
  if not ax:
    f = plt.figure()
    
    ax  = f.add_subplot(111, projection = '3d')
  if scatter:  
        ax.scatter(xs, ys, zs,s = sizes, **kwargs) 
  else:
        ax.plot(xs, ys, zs, **kwargs) 

