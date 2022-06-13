
from math import acos, pi, sqrt
from numpy import array, asfarray, dot, ndarray,asfarray
from Engine import PointMass

from Engine.Connection import Connection
from .Lattice import Lattice
from typing import List, Optional, Tuple


class System():
  def __init__(self,timeMultiplier:float =1,g:float = 9.81):
    self.timeMultiplier:float=timeMultiplier
    self.g:float = g
    self.lattices:List[Lattice]=[]
    self.gravityDirection=asfarray([0,1])

  def addLattice(self,lattice:Lattice):
    lattice.computeMomentOfInertia()
    self.lattices.append(lattice)

  def update(self,deltaT):
    #open frame
    for lattice in self.lattices:
      lattice.startNewFrame()
      if self.g!=0:
        for point in lattice.points.values():
          impulse=self.gravityDirection*point.mass*deltaT*self.timeMultiplier
          lattice.addImpulse(point,impulse)

    #resolve collisions
    for lattice in self.lattices:
      for otherLattice in self.lattices:
        if lattice is otherLattice:
          continue
        #check collision radii
        posDif=lattice.absolutePosition-otherLattice.absolutePosition
        if sqrt(posDif.dot(posDif))<sqrt(lattice.collisionRadiusSquared)+sqrt(otherLattice.collisionRadiusSquared):
          colliding,collisionPoints=self.collisionCheck(lattice,otherLattice)
          if colliding and collisionPoints!=None:
            self.resolveCollision(deltaT,collisionPoints[0],lattice,collisionPoints[1],otherLattice)

    #TODO calculate fracturing
      
    #close frame
    for lattice in self.lattices:
      lattice.resolveFrame(deltaT*self.timeMultiplier)

  
  #collisions defined as the crossing of two connections
  def collisionCheck(self,a:Lattice,b:Lattice)->Tuple[bool,Optional[Tuple[Connection,Connection]]]:
    
    for connection in a.getSurfaceConnections():
      for otherConnection in b.getSurfaceConnections():
        if self.doConnectionsIntersect(connection,a,otherConnection,b):      
          return True,(connection,otherConnection)
    return False,None

  @staticmethod
  def doConnectionsIntersect(connectionA:Connection,latticeA:Lattice,connectionB:Connection,latticeB:Lattice)->bool:

    def intersecting(pa:ndarray,pb:ndarray,pc:ndarray,pd:ndarray)->bool:
      check=(pa-pc)*(pb-pd)
      return check[0]<=0 and check[1]<=0

    #get absolute positions
    pointA=latticeA.getAbsolutePosition(connectionA.p1.relativePos) 
    pointB=latticeA.getAbsolutePosition(connectionA.p2.relativePos)
    pointC=latticeB.getAbsolutePosition(connectionB.p1.relativePos)
    pointD=latticeB.getAbsolutePosition(connectionB.p2.relativePos)

    return intersecting(pointA,pointB,pointC,pointD)


  def resolveCollision(self,deltaT:float,connectionA:Connection,latticeA:Lattice,connectionB:Connection,latticeB:Lattice):
    print("collinding",connectionA.ID,connectionB.ID)

    def distanceBetween(pointC:PointMass,latticeC,connectionD:Connection,latticeD:Lattice)->float:
      directionVector=latticeD.getAbsolutePosition(connectionD.p1.relativePos)-latticeD.getAbsolutePosition(connectionD.p2.relativePos)
      positionVector=latticeD.getAbsolutePosition(connectionD.p1.relativePos)
      pointVector=latticeC.getAbsolutePosition(pointC.relativePos)
      numerator=((asfarray([1,-1])*(directionVector)).dot(pointVector))-directionVector[0]*positionVector[1]+directionVector[1]*positionVector[0]
      perpendicularDistance=numerator/(sqrt(directionVector.dot(directionVector)))

      u=directionVector
      v=pointVector-positionVector
      angleA=(acos(u.dot(v)/sqrt(u.dot(u)*v.dot(v))))
      if abs(angleA)>pi/2:
        #distance between connectionD.p1 and pointC
        distance=positionVector-pointVector
        if angleA<0:
          return -sqrt(distance.dot(distance))
        return sqrt(distance.dot(distance))
      v=v+directionVector#vector between pointC and connectionD.p2
      angleB=(acos(u.dot(v)/sqrt(u.dot(u)*v.dot(v))))
      if abs(angleB)>pi/2:
        distance=positionVector-pointVector-directionVector
        if angleB<0:
          return -sqrt(distance.dot(distance))
        return sqrt(distance.dot(distance))

      return perpendicularDistance

    #step back [DISABLED]
      collisionResolutionStepSize=0.1#fraction of frame that sould be stepped backkj 
      backTimeStep=-collisionResolutionStepSize*deltaT*self.timeMultiplier
      while self.doConnectionsIntersect(connectionA,latticeA,connectionB,latticeB):
        #resolveNetValues is reversible with negative time
        latticeA.resolveNetValues(backTimeStep)
        latticeB.resolveNetValues(backTimeStep)
    
    #find the collision normal
    dap1=distanceBetween(connectionA.p1,latticeA,connectionB,latticeB)
    dap2=distanceBetween(connectionA.p2,latticeA,connectionB,latticeB)
    dbp1=distanceBetween(connectionB.p1,latticeB,connectionA,latticeA)
    dbp2=distanceBetween(connectionB.p2,latticeB,connectionA,latticeA)

    #normal is a unit vector
    if min(abs(dap1),abs(dap2))<min(abs(dbp1),abs(dbp2)):
      #normal is perpendicular to connection B
      absoluteDirectionVector=array(dot(latticeB.rotationMatrix,connectionB.getHeading()))
      connectionAtA=False
    else:
      #normal is perpendicular to connection A
      #multiply by -1 to make normal point to A
      absoluteDirectionVector=-1*array(dot(latticeA.rotationMatrix,connectionA.getHeading()))
      connectionAtA=True

    normal=asfarray([absoluteDirectionVector[1],-1*absoluteDirectionVector[0]])*1/(sqrt(absoluteDirectionVector.dot(absoluteDirectionVector)))
    #resolve perfectly inelastic collision
    parallelMomentum:float=latticeA.velocity.dot(normal)*latticeA.mass+latticeB.velocity.dot(normal)*latticeB.mass
    newParallelVelocity:float=parallelMomentum*(1/(latticeA.mass+latticeB.mass))
    parallelImpulseForA=latticeA.mass*(newParallelVelocity-latticeA.velocity.dot(normal))
    impulseForA=parallelImpulseForA*normal
    impulseForB=-1*impulseForA
    if connectionAtA:
      latticeA.addImpulse(connectionA.p1,impulseForA/2)
      latticeA.addImpulse(connectionA.p2,impulseForA/2)
      if dbp1<dbp2:
        latticeB.addImpulse(connectionB.p1,impulseForB)
      else:
        latticeB.addImpulse(connectionB.p2,impulseForB)
    else:
      latticeB.addImpulse(connectionB.p1,impulseForB/2)
      latticeB.addImpulse(connectionB.p2,impulseForB/2)
      if dap1<dap2:
        latticeA.addImpulse(connectionA.p1,impulseForA)
      else:
        latticeA.addImpulse(connectionA.p2,impulseForA)


