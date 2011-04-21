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
import compbio.utils.colors as mycolors
import compbio.utils.seismic as seismic

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
    description_cols = split_re.split(desc_open.readline().strip()) + ['Exp_Index']
    description_vals = [split_re.split(l.strip()) for l in desc_open.readlines()]
    for idx, d in enumerate(description_vals): d.append(idx)
    

    data_open = open(data_path)
    weight, tf, exp = zip(*[array(split_re.split(l.strip()), float) 
                           for l in data_open.readlines()])
    exp  = [ e -1 for e in exp]
    description = {}
    for i in range(len(description_cols)): 
      description[description_cols[i]] = [d[i] for d in description_vals]
      
    
    ntf = np.max(tf) + 1
    nexp = len(description.values()[0]) 
    
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
    general = lambda x: x['Perturbations'] == 'NA' and x['Time']  == 'NA' \
      and x['OverexpressedGenes'] == 'NA' and x['DeletedGenes'] == 'NA',
    general_ts = lambda x: x['Perturbations'] == 'NA' and x['Time']  != 'NA' \
      and x['OverexpressedGenes'] == 'NA' and x['DeletedGenes'] == 'NA',
    drug = lambda x: x['Perturbations'] != 'NA' and x['Time']  == 'NA' \
      and x['OverexpressedGenes'] == 'NA' and x['DeletedGenes'] == 'NA',
    drug_ts = lambda x: x['Perturbations'] != 'NA' and x['Time']  != 'NA' \
      and x['OverexpressedGenes'] == 'NA' and x['DeletedGenes'] == 'NA',
    genetic = lambda x: x['Perturbations'] == 'NA' and x['Time']  == 'NA' \
      and ( x['OverexpressedGenes'] != 'NA' or x['DeletedGenes'] != 'NA'),
    genetic_ts = lambda x: x['Perturbations'] == 'NA' and x['Time']  != 'NA' \
      and ( x['OverexpressedGenes'] != 'NA' or x['DeletedGenes'] != 'NA'),
    drug_genetic = lambda x: x['Perturbations'] != 'NA' and x['Time']  == 'NA' \
      and ( x['OverexpressedGenes'] != 'NA' or x['DeletedGenes'] != 'NA'),
    drug_genetic_ts = lambda x: x['Perturbations'] != 'NA' and x['Time']  != 'NA' \
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


def dsi_boxplot(num = 1 ,  method = 'tree', reset = False,
                plot_kcs = True,
                bp_means = False,
                bp_zeros = True, zero_ofs = 1e-6,
                bp_logs = True,
                show_kos = True,
                filter_rows_and_cols = True):

  grid, descriptions = parseNet(num= num, method = method, reset = reset)
  grid = array(grid)
  descriptions = dict(descriptions)
  new_descriptions = {}

  if filter_rows_and_cols:
    #Filter out bad rows and columns
    good_exps = nonzero(np.max(grid,0))[0]

    tf_new_idxs = list(argsort(np.max(grid,1))[::-1])
    new_grid = grid[tf_new_idxs]
    good_tfs = nonzero(np.max(new_grid,1))[0]

    
    #Relabel the descriptions to take filtration into account
    #Assumed that one based indexing may be causing havoc so subtract one from the group.
    for k, value in descriptions.iteritems():
      if 'Genes' in k:
        new_descriptions[k] = [re.sub(re.compile('(\d+)'),\
                                        lambda x:  int(x.group()) in tf_new_idxs and str(tf_new_idxs.index(int(x.group()))) or x.group(), g) 
                               for g in value]
      else:
        new_descriptions[k] = value
      new_descriptions[k] = list(array(new_descriptions[k])[good_exps])
      
    new_grid = new_grid[good_tfs, :]
    new_grid = new_grid[ :,good_exps]
    
    grid = new_grid
    descriptions = new_descriptions


  #Make lambdas to split experiments into categories
  col_choosers = sg_choosers()
  #Split experiments
  exps = {}
  for k, v in col_choosers.iteritems():
    vs = [ dict(zip(descriptions.keys() , elt))
          for elt in  zip(*descriptions.values()) ]    
    exps[k] = nonzero( [v(e) for e in vs ])[0]

  '''Remove 'general' as the values wind up being all zeros.'''
  exps.pop('general')
  
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
  
  do_final_bps = True
  kn_exps = {}

  split_ko_ts = False
  
  kn_exps['ko'] = []
  

  
  def getBPS(**kwargs):
    xlabels = []
    nz_frac_std  = []
    nz_frac_mean = []
    nz_val_std   = []
    nz_val_mean  = []
    
    nz_colvals = []

    for k, ecols in exps.iteritems():
      these_knockouts = array([c for c in knockout_cells if c[1] in ecols])
      exp_cells = array([(i,j) for j in ecols for i in arange(shape(grid)[0])])
      if these_knockouts != []:
        kns_found = [c for c in exp_cells 
                     if  np.sum(greater( np.product(c==these_knockouts,1),0),0)]
        kn_exps['ko'] += kns_found

        nokns_found = [c for c in exp_cells 
                       if not np.sum(greater( np.product(c==these_knockouts,1),0),0)]
      else:
        nokns_found = exp_cells

      cexp = [grid[zip(*exp_cells[\
              nonzero(equal(exp_cells[:,1],col))[0]])] \
                         for col in ecols] 
      
      colwise_fracs = [mean(1.*greater(col,0)) for col in cexp]
      colwise_exprs = [mean(col[nonzero(greater(col,0))]) for col in cexp]
      colwise_exprs = [c if not isnan(c) else 0 for c in colwise_exprs]

      nz_colvals.append(colwise_exprs)
      #if k == 'general_ts': raise Exception()

      nz_frac_std.append(std(colwise_fracs)/sqrt(len(colwise_fracs)))
      nz_frac_mean.append(mean(colwise_fracs))
      nz_val_std.append(std(colwise_exprs)/sqrt(len(colwise_exprs)))
      nz_val_mean.append(mean(colwise_exprs))
      
      if isnan(nz_val_mean[-1]): raise Exception()
      
      xlabels.append(k)

    for k, ecells in kn_exps.iteritems():
      ecells = array(ecells)
      nz_frac_std.append(0)
      nz_val_std.append(0)
      nz_frac_mean.append(mean(greater(grid[zip(*ecells)],0)))
      nz_val_mean.append(mean(grid[zip(*ecells[greater(grid[zip(*ecells)],0)])]))
      nz_colvals.append(grid[zip(*ecells[greater(grid[zip(*ecells)],0)])])
      xlabels.append(k)
      
    return xlabels, array(nz_frac_std),array(nz_val_std),array(nz_frac_mean), array(nz_val_mean), [array(cv) for cv in nz_colvals]
  xlabels, nz_frac_std,nz_val_std,nz_frac_mean, nz_val_mean, nz_colvals = mem.getOrSet(getBPS,on_fail = 'compute', reset = reset)
  
  args = [xlabels.index(x) for x in 
          ['general_ts', 'drug', 'drug_ts', 
           'genetic', 'genetic_ts', 'drug_genetic', 'drug_genetic_ts', 'ko']]
  xlabels, nz_frac_std,nz_cal_std,nz_frac_mean,nz_val_mean =\
      array(xlabels)[args],nz_frac_std[args],nz_val_std[args],nz_frac_mean[args],nz_val_mean[args]
  nz_colvals = [nz_colvals[a] for a in args]

  f = plt.figure(0)
  f.clear()


  plot_type = 'dsi_final'
  if plot_type == 'dsi_final':
    margin = .05
    wid0 = .75
    
    ax0 = f.add_axes([margin,margin, wid0 , 1. - 2* margin], title =  'Experminent mean significances: blue (red) lines denote quartiles (media).')
    ax0.boxplot(nz_colvals[0:-1], widths = [.5] * (len(nz_colvals )-1))
    ax0.set_yscale('log')
    ax0.set_xticklabels(xlabels[:-1])
    
    ax1 = f.add_axes([2*margin +wid0, margin, (1 - margin) - (2 * margin + wid0), 1- 2* margin],sharey = ax0, title = 'TF knockout/OE')
    ax1.boxplot(nz_colvals[-1:],widths = .5)
    ax1.set_xticklabels(xlabels[-1:])
  
    f.savefig(config.dataPath('daniel/figs/final_bp_net{0}_{1}.tiff'.format(num, method)))
  
    return
  elif plot_type == 'twoplots':
    nkeys = len(xlabels)
    if show_kos: xi = arange(nkeys)
    else: xi = arange(nkeys -1)
    
    y1 = nz_val_mean[xi]
    s1 =  nz_val_std[xi]
    y2 = nz_frac_mean[xi]
    s2 =  nz_frac_std[xi]
    
    a1 = f.add_subplot(211, ylim =[0, max(y1)+max(s1)], title = 'mean value of nonzero influences\n standard error across experiments')
    a2 = f.add_subplot(212, ylim =[0,max(y2)+ max(s2)], title = 'mean values of fraction nonzero influences\n standard error across experiments' )
    
    colors = mycolors.getct(nkeys)
    wofs = .15
    b1 = a1.bar(xi+wofs,y1,1.-wofs*2, linewidth = 3,color = colors,  ecolor = 'black')
    b2 = a2.bar(xi+wofs,y2,1.-wofs*2, linewidth = 3,color = colors,  ecolor = 'black' )
    p1,c1,b1 = a1.errorbar(xi+.5, y1, yerr = s1,capsize = 15, elinewidth = 4, color = 'black',linewidth = 0, ecolor = 'black')
    p2,c2,b2 = a2.errorbar(xi+.5, y2, yerr = s2,capsize = 15, elinewidth = 4, color = 'black',linewidth =0, ecolor = 'black')
    for c in c1:c.set_alpha(1.)
    for c in c2:c.set_color('black')
    for c in a2.get_children() + a1.get_children():
        try: 
          if not c in [p1,p2]: c.set_linewidth(4)
        except: pass
        continue
    a2.set_xticklabels([])
    for i in xi:
      a2.text( float(i) + .5,0,xlabels[i] , rotation = '-15',size = '16', ha = 'left',va='top')
    f.savefig(config.dataPath('daniel/figs/latest/{1:03d}_{0}.tiff'.format('no_kos' if not show_kos else 'kos', num )),format = 'tiff')
             
  return

def sig_grid(num = 1 ,  method = 'tree', reset = False,
             plot_kcs = True,
             bp_means = False,
             bp_zeros = True, zero_ofs = 1e-6,
             bp_logs = True,
             show_kos = False,
             filter_rows_and_cols = False):



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
  if bp_logs: all_bps = [log(b + zero_ofs) for b in all_bps]
  bp_lzero = log(zero_ofs)

  boxplots = ax3.boxplot([bp for bp in all_bps], widths= .5)
  for p in boxplots.values():
      for e in p: e.set_linewidth(4)    

  #Annotate the boxplot figure
  ann_str = ''
  for i in range(8):
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

    
  plam = lambda: filter_rows_and_cols and 'nonzero_exps_and_tfs_cells_log/'\
      or bp_zeros and not bp_logs and bp_means and 'zeros_means_nolog/'\
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
  
  
  

  mean_xvals = [ mean(all_bps[i][nonzero(greater(all_bps[i],bp_lzero))]) for i in range(len(all_bps))]
  pdfs, xvals = zip(*[histogram(x, bins=50, range=[-15,8], normed=False) for x in all_bps])
  import compbio.utils.colors as colors
  c = colors.getct(len(pdfs))
  f3 = plt.figure(3)
  f3.clear()
                                   
  sax = f3.add_subplot('111')
  seismic.seismic([array(x,float)/ sum(x) for x in pdfs], xax = xvals[0][:-1],stacked = False, colors = c, xmarkpts = mean_xvals, ax = sax)
  

  f4 = plt.figure(4)
  f4.clear()
  ax = f4.add_subplot('121')
  ax.set_title('(log base 10) of Percentage Nonzero for Experiment Classes')
  percs = log10(array([100*float(len(nonzero(greater(x,bp_lzero))[0])) / len(x) for x in all_bps]))
  ax.plot(percs,linewidth = 6)
  ax.set_yticks(percs)

  names = exps.keys() + ['TF Knockout/OE']
  ax.set_yticklabels(['{1}\n{0}'.format('%2.2f' % (10**p), names[idx]) for idx,p in enumerate(percs)])
  
  ax2 = f4.add_subplot('122')
  ax2.set_title('Mean of Nonzero Experiments for Experiment Classes')
  means = array([mean(bp[nonzero(greater(bp,bp_lzero))]) for bp in all_bps])
  ax2.plot(arange(1,9), means,linewidth = 6)
  ax2.boxplot( [bp[nonzero(greater(bp,bp_lzero))] for bp in all_bps], widths = .5)
  ax2.set_yticks(means)

  names = exps.keys() + ['TF Knockout/OE']
  ax2.set_yticklabels(['{1}\n{0}'.format('%2.2f' % (p), names[idx]) for idx,p in enumerate(means)])
  
  

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
  
  
  
