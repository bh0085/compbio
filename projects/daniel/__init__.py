import compbio.config as config
import re, os
from numpy import *
import numpy as np
import compbio.utils.heatmap as hm
import numpy as np
import compbio.utils.memo as mem
import itertools as it
import matplotlib.pyplot as plt
import scipy.signal as ss
import matplotlib


def figs():
  print 'Drawing figs'
  [sig_grid(method = method, num =num, reset = True, bp_means = True, bp_zeros = True, bp_logs = False) for method in ['tree', 'svm'] for num in [1,3]]
  print '1 done'
  [sig_grid(method = method, num =num, reset = True, bp_means = True, bp_zeros = False, bp_logs = False) for method in ['tree', 'svm'] for num in [1,3]]

  print '2 done'
  [sig_grid(method = method, num =num, reset = True, bp_means = True, bp_zeros = False, bp_logs = True) for method in ['tree', 'svm'] for num in [1,3]]

  print '3 done'
  [sig_grid(method = method, num =num, reset = True, bp_means = False, bp_zeros = True, bp_logs = False) for method in ['tree', 'svm'] for num in [1,3]]

  print '4 done'
  [sig_grid(method = method, num =num, reset = True, bp_means = False, bp_zeros = False, bp_logs = False) for method in ['tree', 'svm'] for num in [1,3]]
  
  print '5 done'
  [sig_grid(method = method, num =num, reset = True, bp_means = False, bp_zeros = False, bp_logs = True) for method in ['tree', 'svm'] for num in [1,3]]
  print 'done'

  
def parseNet(num =  1,method = 'tree', reset = False):
  '''
  Get one of daniel's nets. Allowable numbers are 1-3 and allowable
  types are 'tree', 'svm'
'''

  def setNet(**kwargs):
    method =kwargs.get('method', 'tree')
    num = kwargs.get('num', 1)

    description_path = config.dataPath('::daniel/net%s_chip_features.tsv') % num
    data_path = config.dataPath('::daniel/informativeness/%s%s.txt') %(method,num)
    split_re = re.compile('\s')
    
    desc_open = open(description_path)
    description_cols = split_re.split(desc_open.readline().strip())
    description_vals = [split_re.split(l.strip()) for l in desc_open.readlines()]
    
    data_open = open(data_path)
    weight, tf, exp = zip(*[array(split_re.split(l.strip()), float) 
                           for l in data_open.readlines()])
    description = {}
    for i in range(len(description_cols)): 
      description[description_cols[i]] = [d[i] for d in description_vals]
      
    
    ntf = np.max(tf) + 1
    nexp = len(description.values()[0]) + 1
    
    grid = zeros((ntf,nexp))
    for vals in zip(weight,tf,exp): grid[vals[1], vals[2]] = float(vals[0])
    
    return grid, description
  return mem.getOrSet(setNet,
    reset = reset, 
    register = method,
    name = '%s%s' %(method,num),
    method = method,
    num = num)
  
  
def sg_choosers():
  return dict(
    general_exp = lambda x: x['Perturbations'] == 'NA' and x['Time']  == 'NA' \
      and x['OverexpressedGenes'] == 'NA' and x['DeletedGenes'] == 'NA',
    general_ts = lambda x: x['Perturbations'] == 'NA' and x['Time']  != 'NA' \
      and x['OverexpressedGenes'] == 'NA' and x['DeletedGenes'] == 'NA',
    drug_perturbation = lambda x: x['Perturbations'] != 'NA' and x['Time']  == 'NA' \
      and x['OverexpressedGenes'] == 'NA' and x['DeletedGenes'] == 'NA',
    drug_perturbation_ts = lambda x: x['Perturbations'] != 'NA' and x['Time']  != 'NA' \
      and x['OverexpressedGenes'] == 'NA' and x['DeletedGenes'] == 'NA',
    genetic_perturbation = lambda x: x['Perturbations'] == 'NA' and x['Time']  == 'NA' \
      and ( x['OverexpressedGenes'] != 'NA' or x['DeletedGenes'] != 'NA'),
    genetic_perturbation_ts = lambda x: x['Perturbations'] == 'NA' and x['Time']  != 'NA' \
      and ( x['OverexpressedGenes'] != 'NA' or x['DeletedGenes'] != 'NA'),
    genetic_perturbation_drug= lambda x: x['Perturbations'] != 'NA' and x['Time']  == 'NA' \
      and ( x['OverexpressedGenes'] != 'NA' or x['DeletedGenes'] != 'NA'),
    general_perturbation_drug_ts = lambda x: x['Perturbations'] != 'NA' and x['Time']  != 'NA' \
      and ( x['OverexpressedGenes'] != 'NA' or x['DeletedGenes'] != 'NA' ),
    )
  
def sg_big_hm_annotations(f, ax_box):
  ax0 = f.add_axes([0,0,1,1])
  ax0.annotate( 'Heatmaps of experimental significance in TF-Gene edge prediction',
                ax_box[0:2] + ax_box[2:4], xycoords = 'axes fraction',
                ha = 'right', va = 'bottom',
                xytext = [0,10], textcoords= 'offset pixels',
                size = 'xx-large')
  ax0.annotate( 'TFs',
                ax_box[0:2] + [0,ax_box[3]], xycoords = 'axes fraction',
                ha = 'right', va = 'top',rotation = 90,
                xytext = [-10,0], textcoords= 'offset pixels',
                size = 'x-large')
  ax0.annotate( 'Experiments',
                ax_box[0:2] + [0,ax_box[3]], xycoords = 'axes fraction',
                ha = 'left', va = 'bottom',
                xytext = [0, 10], textcoords= 'offset pixels',
                size = 'x-large')

def sig_grid(num = 1 ,  method = 'tree', reset = False,
             plot_kcs = False,
             bp_means = False,
             bp_zeros = False,
             bp_logs = True):
  grid, descriptions = parseNet(num= num, method = method, reset = reset)

  #Make lambdas to split experiments into categories
  col_choosers = sg_choosers()
  #Split experiments
  exps = {}
  for k, v in col_choosers.iteritems():
    vs = [ dict(zip(descriptions.keys() , elt))
          for elt in  zip(*descriptions.values()) ]    
    exps[k] = nonzero( [v(e) for e in vs ])[0]

  #Make and annotate the heatmap figure
  f = plt.figure(1, facecolor = 'w')
  f.clear()
  axdims= .9
  ax_box = array([.05,.05,axdims,axdims])
  sg_big_hm_annotations(f, ax_box)

  #Set up the sizes of each group axis in the heatmap figure
  kwts = float(sum([len(v) for  v in exps.values()]))
  mwidth = .015
  msize = mwidth*kwts
  kw_total = kwts +  ( msize * (len(exps)-1))
  ofs = 0

  #Mark experiments that knock out TFS
  tf_kn_matches =[ sorted(list(it.chain(\
          nonzero([ 'G{0},'.format(t) in x+',' 
                    for x in  descriptions['DeletedGenes'] ])[0],
          nonzero([ 'G{0},'.format(t) in x+',' 
                    for x in  descriptions['OverexpressedGenes'] ])[0])))
                   for t in range(shape(grid)[0])]
  knockout_tfs = nonzero([len(k) for k in tf_kn_matches])[0]
  knockout_cells = array(list(it.chain(*[ [(i, exp) for exp in tf_kn_matches[i] ] 
                               for i in range(len(tf_kn_matches))])))
  knockout_vals = grid[zip(*knockout_cells)]
  allow_tf_kn = False
  if not allow_tf_kn: grid[zip(*knockout_cells)] = 0

  #Some more heatmap configuration.
  saturation = [np.percentile(grid[nonzero(greater(grid,0))],10),
                np.percentile(grid[nonzero(greater(grid,0))],90)]
  tf_srt = argsort(np.mean(grid,1))
  all_bps = []
  expsums = [np.mean( grid.T[v,:], 1) for v in exps.values()]
  max_sum = np.max((list(it.chain(*expsums))))

  #For each experiment class, plot a heatmap and overlay per exp sums
  for k , v in exps.iteritems():
    #Axes positioning
    wid = len(v)
    ax_ofs =  array([ofs/kw_total, 0, (wid) / kw_total,1.])
    ax_box = array([.05,.05,0.,0.])
    ax_ofs = (ax_ofs * axdims) + ax_box

    #Make heatmap axes.
    ax = f.add_axes(ax_ofs, frameon = False)
    sums = np.mean(grid.T[v,:],1)
    exp_srt = argsort(sums)[::-1]
    hm.heatMap( grid.T[v[exp_srt],:][:,tf_srt], axes = ax,
                vmin = saturation[0],
                vmax = saturation[1])

    #Make overlay axes.
    ax2 = f.add_axes(array(ax_ofs) +  array([0,0,0,0]),
                     frameon = True,
                     axisbg = 'none',
                     xticks = [],
                     yticks = [])
    
    #Make the axes look the way I like em
    for a in ax2.spines.values():
      a.set_linewidth(2)
      a.set_alpha(.5)
    these_knockouts = nonzero([c [1]in v for c in knockout_cells])
    kc = knockout_cells[these_knockouts]
    kv = knockout_vals[these_knockouts]
    
    #If plot kcs is selected, plot the cells corresponding to TF deletion/OE
    if plot_kcs:  
      if len(kc) > 0:
        ax.scatter(*zip(*[( list(v).index(x[1]),x[0]) for x in kc]), s =50, 
                  color = 'none', edgecolor = 'black', linewidth = 3)
    color = 'blue'
    ax2.plot(sums[exp_srt],
            linewidth = 4, color = color)

    if bp_means: bpelts = sums
    else: bpelts = grid.T[v,:].flatten()
    if not( bp_zeros ): bpelts = bpelts[nonzero(bpelts)]
    all_bps.append(bpelts)

    ax2.set_xlim([0,wid])
    ax2.set_ylim([0,max_sum])
    ax.set_xlim([0,wid])
    ax.set_ylim([0,shape(grid)[0]])
    ax2.set_xticks([])

    #Annotate each axios
    tbb = matplotlib.transforms.Bbox(ax2.bbox).translated(0,-20)
    t = ax2.text(-2,0, k, 
                 va = 'bottom', ha = 'right',
                 rotation = 90, color = 'black',
                 size = 'x-large', family = 'serif')
    ofs +=  wid + msize


  #Make the boxplot figure
  f2 = plt.figure(3)
  plt.clf()
  

  if bp_means:  bp_kos =  array([  mean(grid.T[g[0],:],0) 
                             for g in it.groupby(sorted(\
        [ko[1] for ko in knockout_cells]))
                             ])
  else: bp_kos = array(knockout_vals)
  if not bp_zeros: bp_kos = bp_kos[nonzero(bp_kos)]

  all_bps = all_bps +  [bp_kos]

  ax3 = f2.add_subplot('111')
  if bp_logs: all_bps = [log(b) for b in all_bps]
  boxplots = ax3.boxplot([bp for bp in all_bps], widths= .5)
  for p in boxplots.values():
      for e in p: e.set_linewidth(4)    

  #Annotate the boxplot figure
  ann_str = ''
  for i in range(9):
    ann_str += '{0}: {1}\n'.format(i+1, (exps.keys() + ['TF Knockout/OE'])[i])
  ax3.annotate(ann_str, [0,1],xycoords = 'axes fraction',
               xytext = [10,-10], textcoords = 'offset pixels',
               va = 'top', ha = 'left')
  ax3.set_title('''Boxplot of significances per experiment type for {3} learning method, Net {4} 

Filtered out were {0} cells corresponding to {1} TFs Knocked out or OverExpressed.
{2} of these cells have nonzero importance and are plotted at x=9,

Showing Means: {5}, Showing zeros: {6}, Plotting logs {7}'''.\
                  format(len(knockout_cells), len(knockout_tfs),
                         len(nonzero(knockout_vals)[0]), 
                         method, num,
                         bp_means, bp_zeros, bp_logs))
  ax3.set_ylabel('significance')
  ax3.set_xlabel('experiment class')
  
  f.savefig(config.dataPath('daniel/figs/{0}_net{1}_heatmaps.tiff'.format(method, num)),
            format = 'tiff')

    
  plam = lambda: bp_zeros and not bp_logs and bp_means and 'zeros_means_nolog/'\
      or not bp_zeros and bp_means and not bp_logs and 'nozeros_means_nolog/'\
      or not bp_zeros and bp_means and bp_logs and 'nozeros_means_log/'\
      or bp_zeros and not bp_means and not bp_logs and 'zeros_cells_nolog/'\
      or not bp_zeros and not bp_means and not bp_logs and 'nozeros_cells_nolog/'\
      or not bp_zeros and not bp_means and bp_logs and 'nozeros_cells_log/'

  dataDir = config.dataPath('daniel/figs/{2}{0}_net{1}_boxplots.tiff'.\
                              format(method, num,plam()))
  print 'saving {0}'.format(dataDir)
  if not os.path.isdir(os.path.dirname(dataDir)): os.mkdir(os.path.dirname(dataDir))
  if os.path.isfile(dataDir): os.remove(dataDir)
  f2.savefig(dataDir,    format = 'tiff')
  
  
  
    

def linearize(array):
  return reshape(array, product(shape(array)))

def compareall():
  f = plt.figure(0)
  plt.clf()
  for i in range(3):
    num = [1,3,4][i]
    
    ax = f.add_subplot('23%s'%(i+1))
    d1,corr1 = compare(num, ax = ax)
    ax2 = f.add_subplot('23%s'%(i+4))
    d2,corr2= compare(num,rnd_tfs = True, ax = ax2)
    ax.set_title('Net %s correlation svm vs. tree importance :\n  %s' % (num, corr1))
    ax2.set_title('Net %s (randomized) \n correlation: %s' % (num, corr2))
    ax.set_xlabel('tree')
    ax.set_ylabel('svm')
def compare(num = 1, reset = False, 
            ax = None, 
            rnd_tfs = False, rnd_exps  = False):

  if not ax:
    plt.clf()
    f = plt.figure(0)
    ax = f.add_subplot('111')
    
  
  net_svm, descriptions = parseNet(num = num, method = 'svm', reset = reset)
  net_tree, descriptions  = parseNet(num = num, method = 'tree', reset = reset)

  if shape(net_svm) != shape(net_tree):

    print '''There is something weird about this data,
the shapes do not match'''
    print '...fixing?'
    if shape(net_svm) > shape(net_tree):
      net_tree2 = zeros(shape(net_svm))
      inds = list(it.product(*[arange(s) for s in shape(net_tree)]))
      net_tree2[zip(*inds)] =\
          net_tree[zip(*inds)]
      net_tree = net_tree2
    else:
      net_svm2 = zeros(shape(net_tree))
      inds = list(it.product(*[arange(s) for s in shape(net_svm)]))
      net_svm2[zip(*inds)] =\
          net_svm[zip(*inds)]
      net_svm = roll(net_svm2,0)

  if rnd_tfs:
    net_svm = net_svm[random.permutation(range(len(net_svm)))]
  if rnd_exps:
    net_svm = net_svm[:, random.permutation(range(shape(net_svm)[1]))]

    
  net_svm = linearize(net_svm)
  net_tree = linearize(net_tree)

  srted = argsort(reshape(net_tree,product(shape(net_tree))))
  srted = srted[nonzero(greater(net_svm[srted] * net_tree[srted], 0))[0]]
  


  print 'Correlation coefficient: '
  corr =    corrcoef(log(net_svm[srted]), log(net_tree[srted]))
  print corr
  ax.plot(log(net_svm[srted]), log(net_tree[srted]),
          marker = 'o',
          linestyle = '',
          color = 'red',
          alpha = .2)
  return [[net_svm, net_tree, srted],corr[0,1]]
  #hm.heatMap(grid, xlabel = 'TFs', ylabel = 'Chips')
  
  
  
