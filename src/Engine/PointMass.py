from dataclasses import dataclass
from typing import List
from numpy import ndarray,zeros

@dataclass
class PointMass:
  __slots__=("relativePos","relativeVelocity","mass","connections","ID")
  
  def __init__(self,mass:float,relativePos:ndarray,ID:int):
    self.mass:float=mass
    self.relativePos:ndarray=relativePos
    self.relativeVelocity:ndarray=zeros([2])
    self.ID:int=ID
  
  def __hash__(self) -> int:
    return self.ID

  def __eq__(self, __o: object) -> bool:
    if isinstance(object,PointMass):
      return self.ID==__o.ID  and self.mass==__o.mass and self.relativePos==__o.relativePos and self.relativeVelocity == __o.relativeVelocity# type: ignore
    else:
      return False