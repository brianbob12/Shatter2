from dataclasses import dataclass
from typing import List
from numpy import ndarray,zeros

@dataclass
class PointMass:
  __slots__=("relativePos","relativeVelocity","mass","connections","ID")
  
  def __init__(self,mass:float,relativePos:ndarray,ID:int):
    self.mass=mass
    self.relativePos:ndarray=relativePos
    self.relativeVelocity:ndarray=zeros([2])
    self.ID:int=ID
  
  def __hash__(self) -> int:
    return self.ID