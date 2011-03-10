from numpy import *
def plotPredictions(y_actual, 
                    y_predicted, 
                    ax = None, 
                    color = 'blue',
                    do_sort = True):
  assert ax
  assert len(y_actual)== len(y_predicted)

  plot_idxs = arange(len(y_actual))
  if do_sort: plot_idxs = argsort(y_actual)
  ax.plot(y_actual[plot_idxs],
          color = color,
          linewidth = 1,
          linestyle = '-',
          zorder = 1)
  ax.plot(y_predicted[plot_idxs],
          color = color,
          linewidth = 4,
          linestyle = ':',
          zorder = 0)
  
