# draggable rectangle with the animation blit techniques; see
# http://www.scipy.org/Cookbook/Matplotlib/Animations
import numpy as np
from numpy import *
import matplotlib.pyplot as plt
from compbio.biopython import biohilbert as bh
from compbio.utils import hilbert as h
import matplotlib.patches as patch


keylist = '''
"enter": zoom in
"backspace": zoom out
"[,]": move left, right
";": preview
'''


class HSensitive:
  movestep = .05
  def __init__(self, axes, bert,actor_kwargs):
    self.actor_kwargs = actor_kwargs
    self.axes = axes
    self.figure = axes.figure
    self.canvas = axes.figure.canvas
    self.fig_drawn = 0

    self.press = None
    self.background = None
    self.highlight = None
    self.pressed = None

    self.fg_artists = []
    self.preview_artists = None
    self.preview_angle = 86
    self.preview_angle_step = 4

    self.zoomwin = .005
    self.zoomloc = 0
    self.movestep = .05
    self.actor = bert
    self.preview_auto = 0

  def connect(self):
    self.ciddraw = self.canvas.mpl_connect(
      'draw_event' , self.on_draw )
    self.cidpress = self.canvas.mpl_connect(
      'button_press_event', self.on_press)
    self.cidrelease = self.canvas.mpl_connect(
      'button_release_event', self.on_release)
    self.cidmotion = self.canvas.mpl_connect(
      'motion_notify_event', self.on_motion)    
    self.cidkeyboard = self.canvas.mpl_connect(
      'key_press_event', self.on_keypress)

  def zoomin(self):
    self.actor.setActorScope(self.getWindow())
    self.redraw_all()
    
  def zoomout(self):
    self.actor.setActorScope([-.75,1.75])
    self.redraw_all()

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
  
    self.preview_artists = self.actor.preview_zoom(
      self.axes, self.preview_artists,
      frac_window[0],frac_window[1],
      **{'angle':preview_angles[mod(self.preview_angle,
                                    len(preview_angles))],
         'lvl':preview_lvls[mod(self.preview_angle,
                                len(preview_angles))]})    
    for p in self.preview_artists:
      p.set_animated = True
      #p.set_axes(self.axes)

  def preview_less(self):
    self.preview_angle -=1
    self.preview_zoom()
  def preview_more(self):
    self.preview_angle +=1
    self.preview_zoom()

  def update_preview(self):
    if  self.preview_artists == None:
      return
    for p in self.preview_artists:
      self.axes.draw_artist(p)    
      self.canvas.blit(self.axes.bbox)
  def redraw_all(self):


    if self.fig_drawn:
      self.actor.draw_holding(self.axes)

    self.axes.clear()
    self.axes.set_autoscale_on(False)
    self.axes.set_xlim(self.actor.get_xlim())
    self.axes.set_ylim(self.actor.get_ylim())
  
    
    self.actor.draw_bg(self.axes)
    self.canvas.draw()

    self.actor.draw_fg(self.axes,self.fg_artists, **self.actor_kwargs)
    self.actor.draw_overlay(self.axes)
    
    for a in self.fg_artists:
      a.set_animated(True)
      self.axes.draw_artist(a)
      
    self.canvas.blit(self.axes.bbox)
    #self.canvas.draw()
    #return
    self.background = self.canvas.copy_from_bbox(self.axes.bbox)
    self.fig_drawn = 1

  def update_bg(self):
    if not self.background:    
      self.background = self.canvas.copy_from_bbox(self.axes.bbox)
    #restore the original bg
    self.canvas.restore_region(self.background)   


  def update_all(self):
    self.update_bg()
    self.update_highlight()
    if self.preview_auto:
      self.preview_zoom()
    self.update_preview()

  def on_keypress(self, event):
    if event.key == 'enter':
      self.zoomin()
    if event.key == '\\':
      self.zoomout()
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
    if event.key == 'r':
      self.redraw_all()
    if event.key == 'p':
      self.preview_auto = 1-self.preview_auto
    if event.key == 'escape':
      self.disconnect()
    self.update_all()
  def update_highlight(self):
    #make sure that the highlight object is nonnull
    if not self.highlight:
      self.highlight =self.axes.plot([0,1],[0,1],
                                     linewidth = 4,
                                     color = 'red',
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

  def on_draw(self, event):
    if not self.fig_drawn:
      self.redraw_all()

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
    'disdisconnect all the stored disconnection ids'
    self.canvas.mpl_disconnect(self.ciddraw )
    self.canvas.mpl_disconnect(self.cidrelease )
    self.canvas.mpl_disconnect(self.cidmotion)
    self.canvas.mpl_disconnect(self.cidkeyboard )
    self.canvas.mpl_disconnect(self.cidpress )
    #self.figure.clear()
    plt.close()



def show_actor(actor,
               fig = 1,
               actor_kwargs = {}
               ):
  
  f = plt.figure(fig, figsize = [10,10])
  for call in f.canvas.callbacks.callbacks['key_press_event'].items():
    f.canvas.mpl_disconnect(call[0])

  ax = f.add_axes([0,0,1,1])
  ax.set_autoscale_on(False)
  ax.set_xlim(actor.get_xlim())
  ax.set_ylim(actor.get_ylim())
  hs = HSensitive(ax,actor, actor_kwargs)
  hs.connect()
  plt.show()

