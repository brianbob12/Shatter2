from typing import Dict, List,Tuple

from numpy import asfarray, ndarray,zeros,ones,sqrt,pi,array,dot,matmul,identity,float32
from numpy.linalg import inv
from tensorflow import Tensor
from tensorflow import matmul as tfMatmul
from math import sin,cos
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
    self.incidentImpulses:List[Tuple[PointMass,ndarray,ndarray]]=[]
    self.connectionStresses:Dict[Connection,float]={}

  def getAbsolutePosition(self,relativePosition:ndarray)->ndarray: 
    return(array(dot(self.rotationMatrix,relativePosition)+self.absolutePosition))

  #def computeMomentOfInertia(self):
  #  self.momentOfInertia:float=0
  #  for point in self.points:
  #    self.momentOfInertia+=point.mass*(point.relativePos.dot(point.relativePos))

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

  def startNewFrame(self):
    self.incidentImpulses:List[Tuple[PointMass,ndarray,ndarray]]=[]

  def addImpulse(self,point:PointMass,absoluteVectorOfImpulse:ndarray):
    #rotate the absoluteVector for impulse by -angle
    relativeVector:ndarray=array(dot(inv(self.rotationMatrix),absoluteVectorOfImpulse))
    self.incidentImpulses.append((point,relativeVector,absoluteVectorOfImpulse))

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

  def resolveFrame(self,deltaT:float):
    oldVelocities:Dict[PointMass,ndarray]={}
    for point in self.points.values():
      oldVelocities[point]=self.getAbsoluteVelocity(point)

    self.resolveNetValues(deltaT)

    velocityDifferences:Dict[PointMass,ndarray]={}
    for point in self.points.values():
      velocityDifferences[point]=self.getAbsoluteVelocity(point)-oldVelocities[point]

    self.connectionStresses=self.solveConnectionImpulses(velocityDifferences)

  def solveConnectionImpulses(self,velocityDifferences:Dict[PointMass,ndarray])->Dict[Connection,float]:
    pointConnectionsMatrix:ndarray=zeros(shape=[len(self.points),len(self.connections)])
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

      directionVector:ndarray=connection.getHeading()
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

      #TODO torque

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
      changeInAngularVelocity=changeInAngularVelocity-(point.relativePos[0]*impulseRelative[1]+point.relativePos[1]*impulseRelative[0])/self.momentOfInertia

    #step 2 find new velocity and angular velocity
    self.velocity=self.velocity+changeInVelocity
    self.angularVelocity=self.angularVelocity+changeInAngularVelocity

    #step 3 find new lattice pos,velocity and angle
    #rotate the relative velocity ot find the absolute velocity
    self.absolutePosition=self.absolutePosition+self.velocity*deltaT
    self.rotate(self.angularVelocity*deltaT)