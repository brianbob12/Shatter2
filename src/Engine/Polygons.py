from re import L
from typing import Dict,Tuple,List
from Engine.PointMass import PointMass
from Engine.Lattice import Lattice
from numpy import asfarray
import random

def rect(
  x:float,
  y:float,
  height:float,
  width:float,
  mass:float,
  pointDensity:float,
  connectionRadius:float,
  connectionStrength:float
  )->Tuple[Lattice,List[PointMass]]:


  lat=Lattice(asfarray([x,y]))

  area=height*width
  pointsToAdd=int(area*pointDensity)
  massPerPoint=mass/pointsToAdd

  points:Dict[str,List[PointMass]]={}
  edgePoints=[]

  def getKey(px:float,py:float)->str:
    return str(int(px/connectionRadius))+":"+str(int(py/connectionRadius))

  def addPoint(px:float,py:float)->PointMass:
    point=lat.addPointMass(massPerPoint,asfarray([px,py]))
    if not getKey(px,py) in points.keys():
      points[getKey(px,py)]=[] 
    points[getKey(px,py)].append(point)
    return point

  # add perimeter
  for i in range(int(height*pointDensity)):
    py=-height/2 + i/pointDensity +0.5/pointDensity
    #the +0.5/pointDensity centers the line
    edgePoints.append(addPoint(width/2+x,py+y))
    edgePoints.append(addPoint(-width/2+x,py+y))
    pointsToAdd-=2
  for i in range(int(width*pointDensity)):
    px=-width/2 + i/pointDensity+0.5/pointDensity
    edgePoints.append(addPoint(px+x,height/2+y))
    edgePoints.append(addPoint(px+x,-height/2+y) )
    pointsToAdd-=2


  #fill inside
  for i in range(pointsToAdd):
    px=random.random()*width - width/2
    py=random.random()*height - height/2
    addPoint(px+x,py+y)
    
  # add connections
  requiredSquaredDistance=connectionRadius**2
  for key in points.keys():
    sk=key.split(":")
    rx=int(sk[0])
    ry=int(sk[1])
    neighboringRegions=[]
    for i in range(-1,1,1):
      for j in range(-1,1,1):
        k=str(rx+i)+":"+str(ry+j)
        if k in points.keys():
          neighboringRegions.append(k)

    for point in points[key]:
      for regionKey in neighboringRegions:
        for otherPoint in points[regionKey]:
          if point.ID==otherPoint.ID:
            continue 
          delta=point.relativePos-otherPoint.relativePos
          distanceSquared=delta.dot(delta)
          if distanceSquared<requiredSquaredDistance:
            lat.connect(point,otherPoint,connectionStrength)

  return lat,edgePoints
