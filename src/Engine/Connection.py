
from numpy import ndarray,sqrt
from dataclasses import dataclass
from .PointMass import PointMass

@dataclass
class Connection:
  ID:int
  p1:PointMass
  p2:PointMass
  breakingImpulse:float

  #relative heading
  def getHeading(self)->ndarray:
    out =self.p1.relativePos-self.p2.relativePos
    #sqrt(v.dot(v)) is the fastest way to find the abs of a vector https://stackoverflow.com/questions/9171158/how-do-you-get-the-magnitude-of-a-vector-in-numpy
    out/=sqrt(out.dot(out))
    return out

  def getHeadingFrom(self,p:PointMass)->ndarray:
    if p is self.p2:
      return(self.getHeading())
    elif p is self.p1:
      return(self.getHeading()*-1)
    else:
      raise(Exception())

  def __hash__(self) -> int:
    return self.ID