import matplotlib.pyplot as plt
from matplotlib.patches import Circle as circ

def color_legend(fig,colors, labels,
                 frameon = False,
                 markersize = 40,
                 pos = 3):
    n = len(colors)
    patches = []
    dax = fig.add_axes([0,0,.1,.1])
    dax.set_visible(False)
    for i in range(n):
        c = colors[i]
        l = labels[i]
        patches.append(circ((0,0),1,
                         ec = 'none',
                         fc = c,
                         ))

    fig.legend(patches,
               labels,pos,
               frameon = frameon
               )
