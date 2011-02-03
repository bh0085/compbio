from numpy import *


class lilturtle():
  def __init__(self,angle = 90, step = 1.0):
    self.unit = array([1,0],float)
    self.pos  = array([0,1],float)
    turning_angle = pi / 180 * angle;
    self.turn_mat = array([[cos(turning_angle), -sin(turning_angle)],
                           [sin(turning_angle), cos(turning_angle)]])
    self.default_step = step;
    self.positions = []
    self.positions.append(array(self.pos))
    self.units = []
    self.units.append(array(self.unit))

  def rotate(self, plus):
    if plus:
      self.unit = dot(self.turn_mat, self.unit)
    else:
      self.unit = dot(self.unit, self.turn_mat)

  def forward(self, step = None):
    if step == None:
      step = self.default_step
    self.pos += self.unit * step
    self.units.append(array(self.unit))
    self.positions.append(array(self.pos))
  
  def hilbert(self,level,plus = True):
    if level == 0:
      return

    self.rotate(not plus)
    self.hilbert(level - 1, not plus )
    self.forward()
    self.rotate( plus)
    self.hilbert(level - 1, plus )
    self.forward()
    self.hilbert(level - 1, plus )
    self.rotate(plus)
    self.forward()
    self.hilbert(level - 1, not plus)
    self.rotate(not plus)
    
  def curve(self):
    p = array(self.positions)
    return p
  def getUnits(self):
    return array(self.units)

  def N(self, level):
    tot = 0 
    for i in range(level):
      tot = 3 + 4 * tot 

    return tot

  def inverseN(self, N):
    tot = 0
    count = 0
    while tot < N:
      tot = 3 + 4 * tot
      count += 1 
    return count
      
    
