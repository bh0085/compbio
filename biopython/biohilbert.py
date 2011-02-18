import itertools as it
import compbio.utils.hilbert as h
from numpy import *
import matplotlib.pyplot as plt
import scipy as sp
import numpy as np
import compbio.utils.dh2 as dh2
from   matplotlib.collections import LineCollection as lc
import compbio.biopython.bhhash as bhhash


#MUST FIX ALL REFS TO SELF.START!!!
  
class BioHilbert(dh2.dh2_actor):
  def __init__(self, **kwargs ):
    self.angle = kwargs.get('angle',90)
    self.levels = kwargs.get('lvls',1)
    self.skeleton = None
    self.current_window = array([0.,1.])
    self.select_skeleton = True



  def _fracForCoords(self, coords):
      return self.speed_hash.frac(coords)
  def _setFracRange(self, frac_range):
    '''set the globalbase range by inputting a fraction'''
    global_range = [self.globalFrac(frac_range[0]), 
                    self.globalFrac(frac_range[1])]
    self.frac_range = global_range



  #member initialization methods.
  def _makeHilbert(self):
    self.bert = h.Hilbert(self.levels,self.angle)
  def _makeHash(self,**kwargs):
    self.speed_hash = bhhash.BHHTree()
    self.speed_hash.make(self.element_lists,**kwargs)



  def init(self,elements_lists, 
           score_only = False,
           **kwargs):

    #Elemental setup
    self.element_lists = elements_lists
    self.frac_range = [0,1]


    #Lambdas to set eventual plots.
    self.speedfun = kwargs.get('speedfun', lambda x: 1)
    self.colorfun = kwargs.get('colorfun', lambda x: 'red')
    self.widthfun = kwargs.get('widthfun', lambda x: '10')
    #self.zfun = zfun

    #make members
    self._makeHilbert()
    self._makeHash(**kwargs)             

    #misc settings
    self.ndraw = kwargs.get('ndraw', 100)



  #interfaces with the speedhash
  def getDrawElements(self, 
                  score_only = False, 
                  n = -1,
                  local_range = None,
                  score_cutoff = -1):

    global_fracs = [self.globalFrac(local_range[0]),
                    self.globalFrac(local_range[1])]
    elts, addresses = self.speed_hash.eltsInRange( global_fracs, n,)
    return elts, addresses
  

  def draw_elements(self, collection, 
                    hilbert,
                    local_range = [0,1], 
                    draw_scope = None,
                    n = 100,
                    score_cutoff = -1):

    
    elts, addresses = self.getDrawElements(False, 
                                           n, 
                                           local_range,
                                           score_cutoff = score_cutoff)

    
    v0 = hilbert.vertices(0,1)
    verts,linewidths,colors,alphas = [], [], [], []   
    srtorder = argsort(map(lambda x: x['speed'], elts))[::-1]
    
    for i in srtorder:
      e,address = elts[i], list(addresses[i])
      frange = array([self.localFrac(self.fracForElement(e, address, 'start'),draw_scope),
                      self.localFrac(self.fracForElement(e, address, 'end'),draw_scope)])

      frange[less(frange,0)] = 0
      frange[greater(frange,1)] = 1
      
      assert frange[1] - frange[0] > 0
      
      v = hilbert.vertices(frange[0],frange[1], snap =False)

      verts.append(v)
      colors.append(self.colorfun(e))
      linewidths.append(self.widthfun(e))


    collection.set_verts(verts)
    collection.set_color(colors)
    collection.set_linewidths(linewidths)
    
  def add_element(self, element, collection, hilbert):
    
    pass


  #handlers for the current range

  #length conversions from bp to fractional and back

  def fracForElement(self, e, address,  which_end = 'start'):
    address = list(address)
    if which_end == 'start':
      position = e['start']
    else:
      position = e['start'] + e['length']
    address.insert(0,position)
    return self._fracForCoords(address)

  def globalFrac(self, frac, frac_range = None):
    if frac_range == None:
      frac_range= self.frac_range
    '''stretch the fractional value given to the full scale'''
    return (frac * (frac_range[1] - frac_range[0]) +\
              frac_range[0])

  def localFrac(self, frac, frac_range = None):
    if frac_range == None:
      frac_range = self.frac_range
    return (frac - frac_range[0]) / \
        (frac_range[1] - frac_range[0])

  #Actor methods
  def setActorScope(self, local_frac):
    self._setFracRange(local_frac)


  def draw_fg(self,ax, artists,**kwargs):

    qdraw = kwargs.get('qdraw',False)
    if qdraw:
        if len(artists) == 0:
          artists.append(ax.add_collection(lc([zeros((2,2)),zeros((2,2))])))
        self.draw_elements(artists[0],self.bert,
                           n = self.ndraw, 
                           score_cutoff = 1)
    else:
      
      elts, addresses = self.getDrawElements()
      for i in range(len(elts)) :
        e, address = elts[i], list(addresses[i])
        fracs =[self.fracForElement(e,address, 'start'),
                self.fracForElement(e,address, 'end')]
        verts = self.bert.vertices(fracs[0], fracs[1],
                                   snap =False)
        artists[0] = ax.plot(verts[:,0],verts[:,1],
                             color = self.colorfun(e),
                             alpha = self.alphafun(e),
                             zorder= self.zfun(e),
                             linewidth = self.widthfun(e))


    skeleton = kwargs.get('skeleton',-1)
    if not self.skeleton:
      self.skeleton = h.Hilbert(skeleton, self.angle)

    if skeleton != -1:
      skel = self.skeleton
      if kwargs.get('skeloff', 0):
        u = skel.units[1:] * skel.step
        n = array([-u[:,1].T, u[:,0].T]).T
        v = skel.curve[:-1] + n * kwargs.get('skeloff') + u / 2

      else:
        v = skel.vertices(0,1.)

      splot = ax.plot(v[:,0], v[:,1],
                zorder = -1,
                linewidth = kwargs.get('skelwidth', 1), 
                alpha = kwargs.get('skelalpha',1),
                color = 'black')
      if len(artists) == 1:
        artists.append(splot[0])
      else:
        artists[1] = splot[0]
        
    ax.set_xticks([])
    ax.set_yticks([])

  def verts_for_range(self, fracs):
    v = array(self.skeleton.vertices(fracs[0],fracs[1]))
    return v.T

  def frac_for_point(self, ptxy):
    return float(self.bert.closest(ptxy)) / len(self.bert.curve)
  def get_xlim(self):
    return array([-.25,1.25])
  def get_ylim(self):
    return array([-1.25,.25])



  def preview_zoom(self, 
                   axes, artists,
                   frac_start,frac_end,
                   **kwargs):

    #FIRST GET XY COORDINATES FOR THE PREVIEW ENDS
    if self.select_skeleton:
      bert = self.skeleton
    else:
      bert = self.bert
    ptstart = bert.map(frac_start)
    ptend = bert.map(frac_end)


    import matplotlib.lines as l

    refined = h.Hilbert(kwargs.get('lvl',4),kwargs.get('angle',90),
                        start= ptstart,
                        finish =ptend)
    verts = refined.vertices(0,1)
    
    if artists == None or len(artists) == 0:
      artists = [axes.plot(verts.T, 
                               linewidth =  40, color = 'white', alpha = .9)[0]]
    else:
      artists[0].set_xdata(verts[:,0])
      artists[0].set_ydata(verts[:,1])

    if len(artists) < 2:
      artists.append(axes.add_collection(lc((verts, verts))))
    

    #global_draw_scope:
    draw_scope = [self.globalFrac(frac_start),self.globalFrac(frac_end)]
    self.draw_elements(artists[1],
                       refined,
                       [frac_start, frac_end], #coordinates in the local scope
                       draw_scope  = draw_scope,
                       n = 200, 
                       score_cutoff = 1)

    


    return artists
