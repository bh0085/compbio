from numpy import *
import numpy as np
class BHHNode():
  def __init__(self):
    self.child_nodes = None
    self.xvals = None
    self.yvals = None
    


class BHHTree():
  '''note: right now, I assume that if a fractional query is made,
  it will fall inside of a child_node...
  ... the coordinate system that I am using demands this:
  all coordinates but the deepest are 'elementwise' whereas
  the deepest is 'basewise' and can fall anywhere'''

  #for the bhhtree position and fraction finding algorithms
  #contained, I assume that intervals are all stored in the
  #the order of their starting position. This allows me to
  #use searchsorted when getting elements and reduces runtime
  #to log(n).

  def make(self, data,**kwargs):
    '''make an interpolated hash tree for the data given'''
    self.root = self._make(data, **kwargs)
  
  def _make(self, elts, **kwargs):
    ''' elts should be a list of dictionaries. 
    for now, they can either hav e asingle key, 'children'
    or 3 keys: 'start, len, score' 
    
    in the case where start, len, score are specified, the
    hashing algorithm originally from biohilbert is called
    to compute a score proportional mapping of 
    [min(starts), max(starts) + len] onto the full score
    proportional interval.

    in the case of children, _make will be called for the
    children of each element; each call will return a
    dictionary having 'start, len, speed', they are
    merged into elts2 which is then sorted as above. '''
    
    #raise Exception()
    assert ('children' in elts[0].keys()) or \
        ('start' in elts[0].keys())

    mynode = BHHNode()
    #RIGHT NOW, this is a really simple function that just
    #computes length and offset values for the current by
    #summing the lengths of children.

    if 'children' in elts[0].keys():
      mynode.child_nodes =[]
      elt_dicts = []
      ofs = 0
      for e in elts:
        out = self._make(e['children'])
        mynode.child_nodes.append(out)
        
        child_size = max(out.yvals)
        elt_dicts.append(dict(length =child_size,
                              start = ofs,
                              speed = 1))

        ofs += child_size * kwargs.get('hash_spacer', 1.0)
      elts = elt_dicts
      mynode.lookup = 'elts'
    else:
      mynode.lookup = 'bases'

    elts = sorted(elts, key = lambda x: x['start'])
    xswitches, yswitches = [], []
    marker0 , marker1 = 0,0

    speed_stack, start_end_stack = [], []
    speed_ends = {}
    searchlam = lambda x:x[2]

    mynode.ends = array(map(lambda x: x['start'] + x['length'], elts))
    mynode.starts = array(map(lambda x: x['start'] , elts))
    mynode.speeds = array(map(lambda x: x['speed'] , elts))
    mynode.ids = array(map(lambda x: x.get('id', 'blank'),elts))

    for e in elts:
      
      while  len(speed_stack) > 0 and \
            start_end_stack[-1][1] < e['start'] :
      
        speed = speed_stack.pop()
        start,end = start_end_stack.pop()
        if end > marker1:
          xswitches.append(marker1)
          yswitches.append(speed)
          marker1 = end

      etop = False
      if len(speed_stack) == 0 or \
            e['speed'] > speed_stack[-1]:
        etop = True

      if speed_ends.get(e['speed'], 0) > e['length'] + e['start']:
        continue

      insrt = searchsorted(speed_stack, e['speed'])
      speed_stack.insert(insrt ,e['speed'])
      start_end_stack.insert(insrt, [e['start'],e['start']+e['length']])
      speed_ends[e['speed']] = e['start'] + e['length']
      

      if etop:
        xswitches.append(e['start'])
        yswitches.append(e['speed'])
        marker1 = e['start'] + e['length']
    
    
    while len(speed_stack) > 0:
      speed = speed_stack.pop()
      start, end = start_end_stack.pop()
      if end > marker1:
        xswitches.append(marker1)
        yswitches.append(speed)
        marker1 = end

    xswitches.append(end)
    yswitches.append(0)

    #compute an integral
    cum = zeros(len(xswitches))
    for i in range(len(xswitches)):
      if i == 0:
        cum[i] = 0
      else:
        cum[i] = cum[i-1] + \
            ( xswitches[i]  - xswitches[i - 1] ) * yswitches[i-1]

    keep = nonzero(not_equal(xswitches, roll(xswitches,1)))[0]
    cum = cum[keep]
    xswitches = array(xswitches)[keep]
    
    mynode.xvals = xswitches
    mynode.yvals = cum
    mynode.elts = elts
    return mynode

  def position(self, frac,**kwargs):
    '''
    return a tuple of positions.
    reflecting the xvalues of the selected
    fraction at each level.
    '''
    coordinates = []
    self._position(self.root, frac, coordinates, **kwargs)
    return coordinates

  def _position(self, node, frac, positions, **kwargs):
    if node.lookup == 'bases':
      out = self._baseAtFrac(node,frac, **kwargs)
      positions.append(out)
    elif node.lookup == 'elts': 
      cidx, cfrac = self._childAtFrac(node, frac, **kwargs)
      child_node = node.child_nodes[cidx]
      elt = self._position(child_node,    
                           cfrac,
                           positions, 
                           **kwargs)
      positions.append( cidx )
      
  def _childAtFrac(self, node, frac, **kwargs):
    assert node.lookup == 'elts'
    frac_yscl = frac * max(node.yvals)
    idxr = searchsorted(node.yvals, frac_yscl, 'right') - 1  
    xval = node.xvals[idxr]
    yval = node.yvals[idxr]
    
    child_match = searchsorted(node.starts, xval, 'right') -1
    match_end = node.ends[child_match]
    '''raise an exception because the fraction selected does not overlap
    the nearest node'''
    assert match_end > xval
    cidx = child_match
    y_end = interp(match_end,node.xvals,node.yvals)
    frac =  (frac_yscl - yval)/(y_end - yval)
    return cidx, frac


  def _baseAtFrac(self, node, frac, **kwargs):
    assert node.lookup == 'bases'
    return round(self._xAtFrac(node, frac))

  def _xAtFrac(self, node, frac,**kwargs):
    frac_yscl = frac * node.yvals[-1]
    idxr = searchsorted(node.yvals, frac_yscl, 'right') - 1  
    xval = node.xvals[idxr]
    yval = node.yvals[idxr]
    if idxr == len(node.xvals) -1:
      remainder = 0
    elif frac_yscl == yval:
      remainder = 0 
    else:
      idx_next = idxr+1
      remainder = (frac_yscl - yval)* \
        (node.xvals[idx_next] - node.xvals[idxr])/ \
        (node.yvals[idx_next] - node.yvals[idxr])
    if 'efound' in kwargs.keys():
      kwargs['efound'].append(node.elts[idxr])
    return xval + remainder

  def address(self, frac, **kwargs):
    address = []
    self._address(self.root,frac,address, **kwargs)
    return address

  def _address(self, node, frac,address,  **kwargs):
    if node.lookup == 'elts':
      cidx, cfrac = self._childAtFrac(node, frac, **kwargs)
      child_node = node.child_nodes[cidx]
      elt = self._address(child_node,    
                           cfrac,
                           address, 
                           **kwargs)
      address.append( cidx )
    if node.lookup == 'bases':
      #call the base at frac function but do not append to elt.
      #simply append the elt to the kwargs dict
      base = self._baseAtFrac(node,frac, **kwargs)
      
            

  def frac(self, coordinates):
    '''
    given a list of coordinates returns a number between 0 and 1
    reflecting a fractional distance of the coordinate across 
    the root node.
    '''
    f = self._frac(self.root, list(coordinates))
    return f
  def _frac(self, node, coordinates):
    '''class method: 
    given a list of coordinates & a node, return a number between 0 and 1
    reflecting fractional distance across the node'''
    coord = coordinates.pop()
    assert node.child_nodes != None or len(coordinates) == 0

    if node.lookup == 'elts':
      #dig down a
      elt_before = node.child_nodes[coord]
      cfrac = self._frac(elt_before,coordinates)
      frac_yval = (1-cfrac) * node.starts[coord] + (cfrac)*node.ends[coord]
      yval =  interp(frac_yval, node.xvals, node.yvals)
      frac =  yval / node.yvals[-1]
      assert frac <= 1 and frac >= 0
      return frac

    elif node.lookup =='bases':
      yout = interp(coord,node.xvals,node.yvals)
      yfrac = yout / node.yvals[-1]
      assert yfrac <=1 and yfrac >= 0
      return yfrac
  def eltWithAddress(self, frac, approx = True):
    assert approx
    elts = []
    address = self.address(frac, efound = elts)
    elt = elts[0]
    return elt, address
    

    
  def eltsInRange(self,frac_range, n):
    #THIS IS A LINEAR OPERATION IN MEMORY AND SPEED
    addr = []
    elts_found = []
    addrs_found = []
    elts = self._eltsInRange(self.root, frac_range, n,
                             addr,
                             elts_found, addrs_found)
    return elts_found, addrs_found
    
  def _eltsInRange(self, node,frac_range, n, addr,
                   elts_found,addrs_found,
                   contained = False):
    assert not contained

    x_ends = [self._xAtFrac(node, frac_range[0]),
              self._xAtFrac(node, frac_range[1])]
    
    
    int_tops = array(node.ends)
    int_tops[greater(int_tops,x_ends[1])] = x_ends[1]
    int_bottoms = array(node.starts)
    int_bottoms[less(int_bottoms,x_ends[0])]=x_ends[0]
    interval_scores = (int_tops-int_bottoms)*node.speeds
    interval_scores *= greater(interval_scores, 0)
    interval_sum = np.sum(interval_scores)
    #find inds that will be allowed to return one or
    #more elements
    if node.lookup =='elts':
      thr_inds = nonzero(greater(interval_scores, interval_sum/n))[0]
    else:
      thr_inds = nonzero(greater(interval_scores, 0)*greater(node.speeds,1))[0]
    
    if node.lookup == 'elts':
      final_sum = np.sum(interval_scores[thr_inds])
      
      for i in thr_inds:
        delta = node.ends[i] - node.starts[i]
        ofs0 = np.max([0,(x_ends[0] -node.starts[i])/delta])
        ofs1 = np.min([1,(1 - ((node.ends[i]-x_ends[1])/delta))])
        dstr = int(floor(interval_scores[i] * n / final_sum))
        new_addr = list(addr)
        new_addr.insert(0,i)
        self._eltsInRange(node.child_nodes[i],
                          [ofs0, ofs1],
                          dstr, new_addr,
                          elts_found,
                          addrs_found,
                          contained = contained)

    elif node.lookup == 'bases':
      #asrt=argsort(interval_scores)[:-10:-1]
      asrt = argsort(interval_scores[thr_inds])
      for i in  thr_inds[asrt][:-n:-1]:
        elts_found.append(node.elts[i])
        addrs_found.append(addr)
    
      
    
