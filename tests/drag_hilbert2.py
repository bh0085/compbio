# draggable rectangle with the animation blit techniques; see
# http://www.scipy.org/Cookbook/Matplotlib/Animations
import numpy as np
from numpy import *
import matplotlib.pyplot as plt
from compbio.biopython import biohilbert as bh
from utils import hilbert as h
import matplotlib.patches as patch



keylist = '''
"enter": zoom in
"backspace": zoom out
"[,]": move left, right
";": preview
'''
class HSensitive:
  movestep = .05
  def __init__(self, axes, bert):
    self.axes = axes
    self.figure = axes.figure
    self.canvas = axes.figure.canvas

    self.press = None
    self.background = None
    self.highlight = None
    self.preview_artist = None
    
    self.pressed = None
    self.preview_verts = None
    self.preview_angle = 86
    self.preview_angle_step = 4

    self.zoomwin = .005
    self.zoomloc = 0
    self.movestep = .05

    self.actor = bert
    
    # draw everything but the selected patch and store the pixel buffer
    canvas = self.canvas
    axes = self.axes
    self.canvas.draw()
    #self.background = self.canvas.copy_from_bbox(self.axes.bbox)
    self.calls = 0
    global keylist
    print keylist

  def connect(self):

    print 'connection'
    self.cidpress = self.canvas.mpl_connect(
      'button_press_event', self.on_press)
    self.cidrelease = self.canvas.mpl_connect(
      'button_release_event', self.on_release)
    self.cidmotion = self.canvas.mpl_connect(
      'motion_notify_event', self.on_motion)    
    self.cidkeyboard = self.canvas.mpl_connect(
      'key_press_event', self.on_keypress)

  def zoomin(self):
    print 'zooming in'
    
  def zoomout(self):
    print 'zoooming out'


  def moveright(self):
    self.zoomloc += self.movestep
  def moveleft(self):
    self.zoomloc -= self.movestep
  def moveto(self, frac):
    self.zoomloc = frac
  def set_winsize(self, frac):
    self.zoomwin = frac

  def getWindow(self):
    return [ np.max([0,self.zoomloc - self.zoomwin]),
             np.min([1,self.zoomloc + self.zoomwin])]
  def preview_zoom(self):
    frac_window = self.getWindow()
    preview_angles = [90,88,83,70,30,10,0,-10,-30,-70,-83,-88,-90]
    preview_lvls =   [5,4,4,2,1,1,1,1,1,2,4,4,5]
  
    self.preview_verts = self.actor.preview_zoom(
      frac_window[0],frac_window[1],
      **{'angle':preview_angles[mod(self.preview_angle,
                                    len(preview_angles))],
         'lvl':preview_lvls[mod(self.preview_angle,
                                len(preview_angles))]})    

  def preview_less(self):
    self.preview_angle -=1
    self.preview_zoom()
  def preview_more(self):
    self.preview_angle +=1
    self.preview_zoom()

  def update_preview(self):
    if  self.preview_verts == None:
      return
    
    if not self.preview_artist:
      self.preview_artist =[]
      self.preview_artist.append(self.axes.plot([0,0],[1,1],
                                                animated = True,
                                                linewidth = 20,
                                                color = 'white',
                                                alpha = .9, zorder = 1)[0])
      self.preview_artist.append(self.axes.plot([0,0],[1,1],
                                           animated = True,
                                           linewidth = 2,
                                                color = 'black',
                                           alpha = 1, zorder = 2)[0])

    for p in self.preview_artist:
      p.set_xdata(self.preview_verts[0])
      p.set_ydata(self.preview_verts[1])
      #redraw the animated portion
      self.axes.draw_artist(p)    
  
    # and blit just the redrawn area
    self.canvas.blit(self.axes.bbox)
       

  def update_bg(self):
    if not self.background:    
      self.background = self.canvas.copy_from_bbox(self.axes.bbox)
    #restore the original bg
    self.canvas.restore_region(self.background)   


  def update_all(self):
    self.update_bg()
    self.update_highlight()
    self.update_preview()

  def on_keypress(self, event):
    if event.key == 'enter':
      self.zoomin()
    if event.key == '[':
      self.moveleft()
    if event.key == ']':
      self.moveright()
    if event.key == 'backspace':
      self.zoomout()
    if event.key == ';':
      self.preview_less()
    if event.key == "'":
      self.preview_more()
    
    self.update_all()
  def update_highlight(self):
    #make sure that the highlight object is nonnull
    if not self.highlight:
      self.highlight =self.axes.plot([0,1],[0,1],
                                     linewidth = 5,
                                     color = 'white',
                                     animated = True)[0]  
    v = self.getverts()
    self.highlight.set_xdata(v[0])
    self.highlight.set_ydata(v[1])

    #redraw the animated portion
    self.axes.draw_artist(self.highlight)      
    # and blit just the redrawn area
    self.canvas.blit(self.axes.bbox)

  def on_press(self,event):
    self.update_bg()
    self.press = event
    self.pressed = True

  def getpt(self, event):
    frac =  self.actor.frac_for_point([event.xdata, event.ydata])
    return frac
  def getverts(self):
    
    frac_range = self.getWindow()
    verts = self.actor.verts_for_range(frac_range)
    return verts


  def on_motion(self,event):

    #check to ensure that we are dealing with a valid event
    if not self.background or not event.xdata:
      return
    #what kind of a motion is occuring?
    if not self.pressed:
      frac= self.getpt(event)
      self.moveto(frac)
    else:
      delta =.005 +  np.sum(square(array([self.press.x, self.press.y]) 
                            - array([event.x, event.y])))  / (square(500))
      self.set_winsize(np.min([1.,delta]))

    self.update_bg()
    self.update_highlight()

  def on_release(self,event):
    self.pressed = False

  def disconnect(self):
    'disconnect all the stored connection ids'
    self.patch.figure.canvas.mpl_disconnect(self.cidpress)
    self.patch.figure.canvas.mpl_disconnect(self.cidrelease)
    self.patch.figure.canvas.mpl_disconnect(self.cidmotion)


def run():
  fig = plt.figure(1)
  ax = fig.add_subplot(111)

  ax.set_autoscale_on(False)
  ax.set_xlim([0,1])
  ax.set_ylim([0,-1])

  actor = h_container(6)
  actor.draw(ax)
  
  hs = HSensitive(ax, actor)
  hs.connect()

  plt.show()
def show_actor(actor,
               fig = 1,
               actor_kwargs = {}
               ):
  
  f = plt.figure(fig)
  ax = f.add_axes([0,0,1,1])
  ax.set_autoscale_on(False)
  ax.set_xlim(actor.get_xlim())
  ax.set_ylim(actor.get_ylim())

  actor.draw(ax, **actor_kwargs)
  hs = HSensitive(ax,actor)
  hs.connect()

  plt.show()

class h_container():
  def __init__(self,lvls):
    self.hilbert = h.Hilbert(lvls)

  def draw(self,ax):
    ax.plot(self.hilbert.curve[:,0], self.hilbert.curve[:,1])
  
  def verts_for_range(self, fracs):
    v = array(self.hilbert.vertices(fracs[0],fracs[1]))
    return v.T

  def frac_for_point(self, ptxy):
    return float(self.hilbert.closest(ptxy)) / len(self.hilbert.curve)
  def get_xlim(self):
    return array([0,1])
  def get_ylim(self):
    return array([0,-1])
