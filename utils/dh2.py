import matplotlib.patches as patches
import matplotlib.text as text
class dh2_actor():
  def __init__(self,lvls):
    self.hilbert = h.Hilbert(lvls)

  def draw_fg(self,ax):
    pass
  
  def draw_bg(self,ax):
    xlim = self.get_xlim()
    ylim = self.get_ylim()

    print [xlim[0], ylim[0]],\
        xlim[1] - xlim[0],\
        ylim[1] - ylim[0]
    rect = patches.Rectangle([xlim[0], ylim[0]], 
                             xlim[1] - xlim[0],
                             ylim[1] - ylim[0], 
                             color = 'white')

    ax.add_patch(rect)
  def draw_overlay(self, ax):
    pass

  def verts_for_range(self, fracs):
    v = array(self.hilbert.vertices(fracs[0],fracs[1]))
    return v.T

  def frac_for_point(self, ptxy):
    return float(self.hilbert.closest(ptxy)) / len(self.hilbert.curve)
  def get_xlim(self):
    return array([-.25,1.25])
  def get_ylim(self):
    return array([-.25,-1.25])
  def draw_holding(self, ax):
    annotation = ax.annotate('Slowly generating \n a new view',
                             [.5,-.5], 
                             verticalalignment = 'center',
                             horizontalalignment = 'center',
                             size = 50 )

    rect = ax.add_patch(patches.FancyBboxPatch([0,-1]
                                               ,1,1,
                                               linewidth = 20,
                                               facecolor = 'white',
                                               alpha = .9,
                                               boxstyle = 'round4'))
    
    ax.draw_artist(rect)
    ax.draw_artist(annotation)
    ax.figure.canvas.blit(ax.bbox)
    ax.figure.canvas.draw()
    
    return annotation
