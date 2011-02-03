from compbio.utils import lilturtle
from numpy import *
import numpy as np
import itertools as it


def N( level):
  tot = 0 
  for i in range(level):
    tot = 3 + 4 * tot 
  return tot

def inverseN( N):
  tot = 0
  count = 0
  while tot < N:
    tot = 3 + 4 * tot
    count += 1 
  return count

class Hilbert():
  def __init__(self, lvl, angle = 90, 
               start = array([0.0,0.0]),
               finish = array([1.0,0.0])):


    
    step = 1.0 / pow(2,lvl)
    t = lilturtle.lilturtle(angle, step)
    t.hilbert(lvl)
    self.curve = t.curve()
    c0 = self.curve[0,:]
    self.curve = self.curve - c0
    last = self.curve[-1,:]
    x,y = last

    delta_in = finish - start
    theta_in = arctan(delta_in[1]/delta_in[0])
    if delta_in[0] < 0: theta_in = theta_in + pi
    theta = arctan(y/x)
    if x < 0: theta = pi + theta
    dtheta = theta - theta_in

    norm_ratio =sqrt( (np.sum(square([x,y]))) /(np.sum(np.square(delta_in))))

    mat = [[cos(dtheta), sin(dtheta)],
           [-sin(dtheta),cos(dtheta)]] *1/norm_ratio
    
    self.step = step / norm_ratio
    curve2 = array(mat[newaxis,:,0] * self.curve[:,newaxis,0] + 
                    mat[newaxis,:,1] * self.curve[:,newaxis,1]) + \
                    start
                      
    self.curve = curve2
    self.units =t.getUnits()
  

  def closest(self,ptxy):
    idx = argmin(np.sum(square(self.curve - ptxy),1))

    #print np.sum(square(self.curve - ptxy)),0)
    #print ptxy
    return idx

    
  def map(self, fraction, snap = False):
    '''
function map(self, fraction, snap = False)
  maps a fractional coordinate to the current hilbert curve. if snap is
  enabled, maps to the nearest vertex. Else, maps to a point along the
  segment lying between vertices.
'''
    c = self.curve
    l = shape(c)[0]
    dist = fraction * (l -1)
    n0 = floor(dist)
    n1 = ceil(dist)
    x  = remainder(dist,1)
    if snap:
      x = round(x)
    y  = 1 - x
    out = c[n0,:] * y + c[n1,:] * x
    return out
    
  def vertices(self, start, finish, 
               snap = False,
               snap_expand = True):
    '''
function vertices(self, start, finish, snap = False)
  maps fractional coordinates start, finish to the current hilbert
  curve and returns intermediate vertices. If snap is enabled, omits
  segment-intermediate points and returns the nearest vertices to
  start, finish.
'''

    c = self.curve
    l = shape(c)[0]
    ends, idxs, round_up = [], [] , []

    for fraction in [start, finish]:
      dist = fraction * (l -1)
      n0 = floor(dist)
      n1 = ceil(dist)

      x  = remainder(dist,1)
      y  = 1 - x
      ends.append( c[n0,:] * y + c[n1,:] * x)
      idxs.append( n0 )
      round_up.append(round(x))

    if snap and snap_expand:
      round_up = [0,1]

    if not snap:
      idx_r = [idxs[0] + 1, idxs[1]]
    else:
      idx_r = [idxs[0] + round_up[0], idxs[1] + round_up[1]]
      
    inds = arange(idx_r[0],(idx_r[1] + 1))
    vals = []
    if not snap:
      vals.append(ends[0])
    vals.extend(c[idx_r[0]:idx_r[1]+1,:])
    if not snap:
      vals.append(ends[1])
    return array(vals)


default_speedfun = lambda x: x.type == 'tRNA' and 1000 or x.type == 'CDS' and 3 or 0
def gInterp(records,speedfun = default_speedfun, method = 'max'):
  '''
gInterp(records, speedfun = default, method = 'max')
  gInterp sets parametrization speeds across the genome.
  if it is left with the default speed function, it stretches coding sequences by 5 and tRNAs by 100

'''
  src = it.ifilter(lambda x: x.type == 'source', records).next()
  end = src.location.end.position
  speed_hash = ones(end +1) 
  for r in records:
    score = speedfun(r)
    loc = [r.location.start.position, r.location.end.position]

    slc = speed_hash[loc[0]:loc[1]+1]
    if method == 'max':
      slc[less(slc,score)] = score
    else:
      raise Exception('method not yet implemented')
    

  for i in range(1,len(speed_hash)):
    speed_hash[i] += speed_hash[i-1]
  speed_hash/= speed_hash[-1]
  return speed_hash

def bInterp(records, speedfun, rng, method = 'max'):
  '''
  Takes as input a list of dictionaries, each element having a 'label' key.
  And speedfun, describing the speed value for each possible value of the label key.
  And rng, specifying the range of the sequence in question for which to construct
  the speed hash.
  
  In violation of the the pythonic standard (but in keeping with the genbank standard)
  the hash table will have rng[1] - rng[0] + 1 elements.
'''  
  n = rng[1] - rng[0] + 1
  speed_hash = ones(n)
  for r in records:
    score = speedfun(r)
    loc =array( [r['start'] , r['end']] )
    loc_hash = (loc - rng[0]) + [ 0 , 1 ]
    
    slc = speed_hash[loc_hash[0]:loc_hash[1]]
    if method == 'max':
      slc[less(slc,score)] = score
    else:
      raise Exception('method not yet implemented')
    
  for i in range(1,len(speed_hash)):
    speed_hash[i] += speed_hash[i-1]
  speed_hash/= speed_hash[-1]
  return speed_hash
