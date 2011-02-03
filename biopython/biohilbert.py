import itertools as it
import utils.hilbert as h
from numpy import *


class BioHilbert():
  def __init__(self, angle = 90, lvls = 4 ):
    self.angle = angle
    self.levels = lvls
    self.skeleton = None
    self.current_window = array([0.,1.])
    self.select_skeleton = True

  def initGB(self,gb, 
             score_only = False,
             ndraw = -1,):
    src = it.ifilter( lambda x: x.type == 'source', gb.features).next()
    
    self.start = src.location.start.position
    self.end  = src.location.end.position 


    self.speedfun =  lambda x :\
      x['label'] == 'exon' and 50000. \
      or x['label'] == 'CDS' and 1000. \
      or x['label'] == 'gene' and 50. \
      or 1

    self.colorfun =  lambda x :\
      x['label'] == 'exon' and 'red' \
      or x['label'] == 'CDS' and 'orange' \
      or x['label'] == 'gene' and 'blue' \
      or 'black'

    
    self.widthfun =  lambda x :\
      x['label'] == 'exon' and 30.0 \
      or x['label'] == 'CDS' and 20.0 \
      or x['label'] == 'gene' and 5.0 \
      or 1.0

    
    self.alphafun = lambda x: \
        x['label'] == 'exon' and 1.0 \
        or x['label'] == 'CDS' and .5 \
        or x['label'] == 'gene' and  .5 \
        or .01

    self.zfun = lambda x: \
        x['label'] == 'exon' and 1 \
        or x['label'] == 'CDS' and 0 \
        or -1

    self.bert = h.Hilbert(self.levels,self.angle)
    self.elements = map(
      lambda x: {'label': x.type,
                 'start': x.location.start.position,
                 'end':x.location.end.position},
      gb.features)
    self.ndraw = ndraw

  def makeHash(self):
    return h.bInterp(self.elements, self.speedfun, [self.start, self.end])
  def setHash(self, h):
    self.speed_hash = h

  #Called when a speed hash has already been created.
  #Get a list of elements to draw.
  def getDrawIdxs(self, 
                  score_only = False, 
                  n = -1,
                  base_range = None):

    if n == -1:
      n = self.ndraw
    draw_scores = array(map(lambda x: self.speedfun(x), self.elements))
    if not score_only:
      draw_scores *= map(lambda x: x['end'] - x['start'] + 1, self.elements)
    if base_range:
      draw_scores *= map(lambda x: x['end'] > base_range[0] or x['start'] < base_range[1] and 1. or 0.,elements) 
      
    draw_idxs = argsort(draw_scores)[-n:]
    return draw_idxs
    
    
  def fracForBase(self, base, frac_range = [0.,1.]):
    return self.speed_hash[base] \
        * (frac_range[1] - frac_range[0]) \
        + frac_range[0]
  
    

  def preview_refine_between(self, 
                             start, 
                             finish,
                             **kwargs):
    refined = h.Hilbert(kwargs.get('lvl',4),kwargs.get('angle',90),
                        start = start,
                        finish = finish)
    verts = refined.vertices(0,1)
    return verts


  def draw(self,ax, **kwargs):
    for i in self.getDrawIdxs():
      e = self.elements[i]
      if not e['label'] in ['exon', 'gene']:
        pass  #continue
      fracs =[self.fracForBase(e['start']),
              self.fracForBase(e['end'])]
      verts = self.bert.vertices(fracs[0], fracs[1],
                                 snap =False)
      ax.plot(verts[:,0],verts[:,1],
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

      ax.plot(v[:,0], v[:,1],
              zorder = -1,
              linewidth = kwargs.get('skelwidth', 1), 
              alpha = kwargs.get('skelalpha',1),
              color = 'black')

      ax.set_xticks([])
      ax.set_yticks([])

  def verts_for_range(self, fracs):
    v = array(self.skeleton.vertices(fracs[0],fracs[1]))
    return v.T

  def frac_for_point(self, ptxy):
    return float(self.bert.closest(ptxy)) / len(self.bert.curve)
  def get_xlim(self):
    return array([-.1,1.1])
  def get_ylim(self):
    return array([-1.1,.1])
  def preview_zoom(self, 
                   frac_start,frac_end,
                   **kwargs):
    if self.select_skeleton:
      bert = self.skeleton
    else:
      bert = self.bert
    ptstart = bert.map(frac_start)
    ptend = bert.map(frac_end)
    return array(self.preview_refine_between(ptstart,
                                             ptend, **kwargs)).T
