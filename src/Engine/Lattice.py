from re import I
from typing import Dict, List, Optional,Tuple

from numpy import asfarray, ndarray,zeros,ones,sqrt,pi,array,dot,matmul,identity,float32
from numpy.linalg import inv
from tensorflow import Tensor
from tensorflow import matmul as tfMatmul
from math import acos, atan2, sin,cos
from random import randint

from Engine.SGDengine import SGDengine
from .Connection import Connection

from .PointMass import PointMass

class Lattice:
  def __init__(self,pos:ndarray):
    self.absolutePosition:ndarray=pos#position of center of mass
    self.velocity:ndarray=zeros([2])
    self.angle:float=0#in radians between -Pi and Pi
    self.rotationMatrix:ndarray=array([[1.0,0.0],[0.0,1.0]])#identity

    self.angularVelocity:float=0

    self.mass:float=0
    self.momentOfInertia:float=0

    self.points:Dict[int,PointMass]={}
    self.connections:Dict[int,Connection]={}
    self.pointConnectionMap:Dict[PointMass,List[Connection]]={}
    self.incidentImpulses:List[Tuple[PointMass,ndarray,ndarray]]=[]
    self.connectionStresses:Dict[Connection,float]={}

    self.impulsesLastFrame:List[Tuple[PointMass,ndarray,ndarray]]=[]

    #the maximum distance between the center of mass and a point
    self.collisionRadiusSquared:float=0

    #List of connections on the surface
    self.surfaceRootNode:Optional[PointMass]=None
    self.surfaceConnections:Optional[List[Connection]]=None

  def getAbsolutePosition(self,relativePosition:ndarray)->ndarray: 
    return(array(dot(self.rotationMatrix,relativePosition)+self.absolutePosition))

  def computeMomentOfInertia(self):
    self.momentOfInertia:float=0
    for point in self.points.values():
      self.momentOfInertia+=point.mass*(point.relativePos.dot(point.relativePos))

  #NOTE for now, this won't work if the lattice is rotate at all
  def addPointMass(self,mass:float,absolutePosition:ndarray)->PointMass:
    #recompute center of mass 
    l=len(self.points)
    ID:int=randint(0,int(1E20)) 
    while ID in self.points.keys():
      ID=randint(0,int(1E20))     

    if l==0:
      self.absolutePosition=absolutePosition.copy()
    else: 
      absolutePositionShift:ndarray=(1/(mass+self.mass))*(self.absolutePosition*self.mass+absolutePosition*mass)-self.absolutePosition
      self.absolutePosition=self.absolutePosition+absolutePositionShift

      #update all other points with their new positions
      for point in self.points.values():
        point.relativePos=point.relativePos-absolutePositionShift
        if point.relativePos.dot(point.relativePos)>self.collisionRadiusSquared:
          self.collisionRadiusSquared=point.relativePos.dot(point.relativePos)
          self.surfaceRootNode=point
      
      

    p=PointMass(mass,absolutePosition-self.absolutePosition,ID)
    

    self.points[ID]=p
    self.mass+=p.mass
    
    #This works even though the center of mass is moving
    self.momentOfInertia+=p.mass*(p.relativePos.dot(p.relativePos))

    return p

  def connect(self,p1:PointMass,p2:PointMass,breakingForce:float):
    ID:int=randint(0,int(1E20))
    while ID in self.connections.keys():
      ID=randint(0,int(1E20))

    c=Connection(ID,p1,p2,breakingForce)
    self.connections[ID]=c 
    self.connectionStresses[c]=0.0
  
    # map points to connections
    if not p1 in self.pointConnectionMap.keys():
      self.pointConnectionMap[p1]=[]
    if not p2 in self.pointConnectionMap.keys():
      self.pointConnectionMap[p2]=[]
    self.pointConnectionMap[p1].append(c)
    self.pointConnectionMap[p2].append(c)

  def startNewFrame(self):
    self.impulsesLastFrame:List[Tuple[PointMass,ndarray,ndarray]]=self.incidentImpulses
    self.incidentImpulses:List[Tuple[PointMass,ndarray,ndarray]]=[]

  def addImpulse(self,point:PointMass,absoluteVectorOfImpulse:ndarray):
    #rotate the absoluteVector for impulse by -angle
    relativeVector:ndarray=array(dot(inv(self.rotationMatrix),absoluteVectorOfImpulse))
    self.incidentImpulses.append((point,relativeVector,absoluteVectorOfImpulse))

  def areImpulsesTheSame(self,tolerance:float=1e-0)->bool:
    t2=tolerance**2
    if len(self.impulsesLastFrame)!=len(self.incidentImpulses):
      return False
    if len(self.impulsesLastFrame)==0 and len(self.incidentImpulses)==0:
      return True
    for impulse in self.impulsesLastFrame:
      found=False
      for otherImpulse in self.incidentImpulses:
        #NOTE only the point and relative vector need to match

        #check impulse is acting on the same point
        if not impulse[0] is otherImpulse[0]:
          continue
        difference=impulse[1]-otherImpulse[1]
        #check the relative vector is within tolerance
        if difference.dot(difference)>t2:
          continue
        found=True
        break
      if not found:
        return False
    return True

  def rotate(self,theta:float):
    self.angle+=theta
    if self.angle<=-2*pi:
      self.angle=self.angle+4*pi
    elif self.angle>2*pi:
      self.angle=self.angle-4*pi
    self.rotationMatrix=array([
      [cos(self.angle),-sin(self.angle)],
      [sin(self.angle),cos(self.angle)]
    ])

  def getAbsoluteVelocity(self,point:PointMass)->ndarray:
    #absolute velocity(of the lattice) + absolute tangential velocity(of the point)
    return self.velocity+array(dot(self.rotationMatrix,array(dot(asfarray([[0,1],[-1,0]]),point.relativePos))))*self.angularVelocity

  def resolveFrame(self,deltaT:float,sameImpulses:Optional[bool]=None):
    if sameImpulses==None:
      sameImpulses=self.areImpulsesTheSame()

    if not sameImpulses:
      print("resolving internal impulses")
      oldVelocities:Dict[PointMass,ndarray]={}
      for point in self.points.values():
        oldVelocities[point]=self.getAbsoluteVelocity(point)

    self.resolveNetValues(deltaT)

    if not sameImpulses:
      velocityDifferences:Dict[PointMass,ndarray]={}
      for point in self.points.values():
        velocityDifferences[point]=self.getAbsoluteVelocity(point)-oldVelocities[point]

      self.connectionStresses=self.solveConnectionImpulses(velocityDifferences)

  def solveConnectionImpulses(self,velocityDifferences:Dict[PointMass,ndarray])->Dict[Connection,float]:
    #everything should be in absolute space
    pointConnectionsMatrix:ndarray=zeros(shape=[len(self.points),len(self.connections)])
    #absolute velocity deltas
    pointVelocityDeltas:ndarray=zeros(shape=[len(self.points),2])
    connectionDirectionMatrix:ndarray=zeros(shape=[len(self.connections),2])#x ands y
    netExternalPointImpulsesTimesMass:ndarray=zeros(shape=[len(self.points),2])

    for point,relativeVector,absoluteVector in self.incidentImpulses:
      pointIndex=list(self.points.keys()).index(point.ID)
      netExternalPointImpulsesTimesMass[pointIndex]=relativeVector*point.mass

    for connectionIndex,connection in enumerate(self.connections.values()):

      pointConnectionsMatrix[list(self.points.keys()).index(connection.p1.ID)][connectionIndex]=1/connection.p1.mass
      #indicates the force acts in the opposite direction
      pointConnectionsMatrix[list(self.points.keys()).index(connection.p2.ID)][connectionIndex]=-1/connection.p2.mass

      directionVector:ndarray=array(dot(self.rotationMatrix,connection.getHeading()))
      distance:float=sqrt(directionVector.dot(directionVector))
      connectionDirectionMatrix[connectionIndex]=directionVector/distance
    
    for pointIndex,point in enumerate(self.points.values()):
      pointVelocityDeltas[pointIndex]=velocityDifferences[point]

    connectionImpulses:ndarray=zeros(shape=[len(self.connections),1]) 

    def f(x:Tensor)->Tensor:
      #tension and compression
      #Tensors and ndarrays not compatible in np.matmul
      a=tfMatmul(x,ones([1,2]))
      #elementwise multiplication
      connectionImpulsesWithDirection:Tensor=connectionDirectionMatrix*a

      #shape=[connections,2]
      calculatedVelocities:Tensor=tfMatmul(pointConnectionsMatrix,connectionImpulsesWithDirection) + netExternalPointImpulsesTimesMass
      return pointVelocityDeltas - calculatedVelocities  # type: ignore

    connectionImpulses=SGDengine.solve(f,connectionImpulses,0.1,150,0.1)


    out:Dict[Connection,float]={} 
    for connectionIndex,connection in enumerate(self.connections.values()):
      out[connection]=connectionImpulses[connectionIndex][0]

    return out


  def resolveNetValues(self,deltaT:float):
    #step 1 compute change in velocity and angular velocity
    changeInVelocity:ndarray=zeros([2])
    changeInAngularVelocity:float=0
    for point,impulseRelative,impulseAbsolute in self.incidentImpulses:
      changeInVelocity=changeInVelocity+impulseAbsolute/self.mass

      #change in angular velocity = r/I *impulse
      changeInAngularVelocity=changeInAngularVelocity+(point.relativePos[0]*impulseRelative[1]+point.relativePos[1]*impulseRelative[0])/self.momentOfInertia

    #step 2 find new velocity and angular velocity
    self.velocity=self.velocity+changeInVelocity
    self.angularVelocity=self.angularVelocity+changeInAngularVelocity

    #step 3 find new lattice pos,velocity and angle
    #rotate the relative velocity ot find the absolute velocity
    self.absolutePosition=self.absolutePosition+self.velocity*deltaT
    self.rotate(self.angularVelocity*deltaT)
    
  def computeSurfaceConnections(self)->List[Connection]:
    if self.surfaceRootNode==None:
      #not enough points
      self.surfaceConnections=[]
      return []

    #clockwise turn
    def angleBetween(p1:ndarray,p2:ndarray,p3:ndarray)->float:
      va=p1-p2
      vb=p2-p3
      angle=atan2(va[0],va[1])-atan2(vb[0],vb[1])
      if angle<0:
        angle+=2*pi
      elif angle<2*pi:
        angle%=2*pi
      return angle

    #compute second point on permitter
    bestAngle:float=4
    secondPoint:PointMass=self.surfaceRootNode
    bestConnection=None
    for connection in self.pointConnectionMap[self.surfaceRootNode]:
      if self.surfaceRootNode is connection.p1:
        p=connection.p2
      else:
        p=connection.p1
      #angle with tangent with center
      angle=angleBetween(p.relativePos,self.surfaceRootNode.relativePos,self.absolutePosition)
      if angle<0: angle+=2*pi
      elif angle>2*pi:angle%=2*pi
      if angle<bestAngle:
        bestAngle=angle
        secondPoint=p
        bestConnection=connection
      

    addedNodes=[self.surfaceRootNode,secondPoint]
    if bestConnection != None:
      addedConnections:List[Connection]=[bestConnection]

    def addNext(point:PointMass,lastPoint:PointMass,index:int):
      if index==800:
        #print("exitting")
        return
      nextPoint:Optional[PointMass]=None
      nextConnection:Optional[Connection]=None
      minAngle:float=8#above2pi
      for connection in self.pointConnectionMap[point]:
        if point is connection.p1:
          p=connection.p2
        else:
          p=connection.p1
        if p in addedNodes and (not p is self.surfaceRootNode):
          continue
        angle=angleBetween(p.relativePos,point.relativePos,lastPoint.relativePos)
        if angle<minAngle:
          minAngle=angle
          nextPoint=p
          nextConnection=connection
        elif angle==minAngle:
          #chose the closer point
          if nextPoint!=None:
            d1=point.relativePos-nextPoint.relativePos
            d2=point.relativePos-p.relativePos
            if d2.dot(d2)<d1.dot(d1):
              nextPoint=p
              nextConnection=connection


      if nextPoint == None:
        #backtrack
        if index==2:
          #reached beginning
          return
        #print("backtracking")
        #print(len(self.pointConnectionMap[point]))
        addNext(addedNodes[index-1],addedNodes[index-2],index-1)
      elif nextPoint is self.surfaceRootNode:
        if nextConnection!=None:
          addedConnections.append(nextConnection)
          return
      else:
        addedNodes.append(nextPoint) 
        if nextConnection !=None:
          addedConnections.append(nextConnection)
        addNext(nextPoint,point,index+1)

    addNext(secondPoint,self.surfaceRootNode,2)

    self.surfaceConnections=addedConnections
    return addedConnections

  def getSurfaceConnections(self)->List[Connection]:
    if self.surfaceConnections==None:
      return self.computeSurfaceConnections()
    else:
      return self.surfaceConnections

