from math import cos, pi, sin, sqrt
from re import L
from typing import Dict,Tuple,List
from Engine.PointMass import PointMass
from Engine.Lattice import Lattice
from numpy import asfarray
import random

def getKey(px:float,py:float,connectionRadius:float)->str:
  return str(int(px/connectionRadius))+":"+str(int(py/connectionRadius))

def addPoint(px:float,py:float,pointDict:Dict[str,List[PointMass]],lat:Lattice,massPerPoint:float,connectionRadius:float)->PointMass:
  point=lat.addPointMass(massPerPoint,asfarray([px,py]))
  if not getKey(px,py,connectionRadius) in pointDict.keys():
    pointDict[getKey(px,py,connectionRadius)]=[] 
  pointDict[getKey(px,py,connectionRadius)].append(point)
  return point

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

  

  sqrtpd=pointDensity**0.5

  # add perimeter
  for i in range(int(height*sqrtpd)):
    py=-height/2 + i/sqrtpd +0.5/sqrtpd
    #the +0.5/pointDensity centers the line
    edgePoints.append(addPoint(width/2+x,py+y,points,lat,massPerPoint,connectionRadius))
    edgePoints.append(addPoint(-width/2+x,py+y,points,lat,massPerPoint,connectionRadius))
    pointsToAdd-=2
  for i in range(int(width*sqrtpd)):
    px=-width/2 + i/sqrtpd+0.5/sqrtpd
    edgePoints.append(addPoint(px+x,height/2+y,points,lat,massPerPoint,connectionRadius))
    edgePoints.append(addPoint(px+x,-height/2+y,points,lat,massPerPoint,connectionRadius) )
    pointsToAdd-=2


  #fill inside
  for i in range(pointsToAdd):
    px=random.random()*width - width/2
    py=random.random()*height - height/2
    addPoint(px+x,py+y,points,lat,massPerPoint,connectionRadius)
    
  # add connections
  requiredSquaredDistance=connectionRadius**2
  for key in points.keys():
    sk=key.split(":")
    rx=int(sk[0])
    ry=int(sk[1])
    neighboringRegions=[]
    for i in range(-2,2,1):
      for j in range(-2,2,1):
        k=str(rx+i)+":"+str(ry+j)
        if k in points.keys():
          neighboringRegions.append(k)

    for point in points[key]:
      for regionKey in neighboringRegions:
        for otherPoint in points[regionKey]:
          if point.ID==otherPoint.ID:
            continue 
          #to avoid doubling up connections only connect points with higher IDs as p1
          if point.ID<otherPoint.ID:
            continue
          delta=point.relativePos-otherPoint.relativePos
          distanceSquared=delta.dot(delta)
          if distanceSquared<requiredSquaredDistance:
            lat.connect(point,otherPoint,connectionStrength)

  return lat,edgePoints

def circle(
  x:float,
  y:float,
  r:float,
  mass:float,
  pointDensity:float,
  connectionRadius:float,
  connectionStrength:float
  )->Tuple[Lattice,List[PointMass]]:
  lat=Lattice(asfarray([x,y]))

  area=pi*r**2
  pointsToAdd=int(area*pointDensity)
  massPerPoint=mass/pointsToAdd

  points:Dict[str,List[PointMass]]={}
  edgePoints=[]

  sqrtpd=pointDensity**0.5

  # add perimeter
  for i in range(int(2*pi*r*sqrtpd)):
    theta=i/(r*sqrtpd)
    edgePoints.append(addPoint(r*cos(theta)+x,r*sin(theta)+y,points,lat,massPerPoint,connectionRadius))
  pointsToAdd-=int(2*pi*r*sqrtpd)

  #fill inside
  for i in range(pointsToAdd):
    #chose distance weighted by circumference
    distance=r*sqrt(random.random())
    theta=random.random()*2*pi
    addPoint(distance*cos(theta)+x,distance*sin(theta)+y,points,lat,massPerPoint,connectionRadius)

  # add connections
  requiredSquaredDistance=connectionRadius**2
  for key in points.keys():
    sk=key.split(":")
    rx=int(sk[0])
    ry=int(sk[1])
    neighboringRegions=[]
    for i in range(-1,2,1):
      for j in range(-2,2,1):
        k=str(rx+i)+":"+str(ry+j)
        if k in points.keys():
          neighboringRegions.append(k)

    for point in points[key]:
      for regionKey in neighboringRegions:
        for otherPoint in points[regionKey]:
          if point.ID==otherPoint.ID:
            continue 
          #to avoid doubling up connections only connect points with higher IDs as p1
          if point.ID<otherPoint.ID:
            continue
          delta=point.relativePos-otherPoint.relativePos
          distanceSquared=delta.dot(delta)
          if distanceSquared<requiredSquaredDistance:
            lat.connect(point,otherPoint,connectionStrength)

  return lat,edgePoints




  
