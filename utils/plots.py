import matplotlib.pyplot as plt
from matplotlib.patches import Circle as circ

def color_legend(fig,colors, labels,
                 frameon = False,
                 markersize = 40,
                 pos = 0,
                 ax = None):
    n = len(colors)
    patches = []

    if not ax:
        dax = fig.add_axes([0,0,.1,.1])
        dax.set_visible(False)
        trg = fig
    else:
        trg = ax
    for i in range(n):
        c = colors[i]
        l = labels[i]
        patches.append(circ((0,0),1,
                         ec = 'none',
                         fc = c,
                         ))
    
    trg.legend(patches,
               labels,pos,
               frameon = frameon,prop = {'family':'serif'},
               fancybox = True
               )


def maketitle(ax, title, subtitle = None):
    ax.annotate(title,[.05,.9],
                xycoords = 'axes fraction',
                verticalalignment = 'bottom',
                size = 'x-large'
                )
    if subtitle:
        ax.annotate(subtitle, [.05,.9],
                    size = 'large',
                    xycoords = 'axes fraction',
                    xytext = [0,-5],
                    textcoords = 'offset points'
                    ,verticalalignment = 'top')

def label_lr(ax, label):
    ax.annotate(label,[.95,.05],
                verticalalignment = 'bottom',
                horizontalalignment = 'right',
                size = 'large')
