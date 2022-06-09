
import time
import numpy as np
import Renderer
import Engine

myRenderer=Renderer.Renderer(1024,1024)

myLattice=Engine.Lattice(np.asfarray([0,0]))
p1=myLattice.addPointMass(1,np.asfarray([100,-100]))
p2=myLattice.addPointMass(1,np.asfarray([100,100]))
p3=myLattice.addPointMass(1,np.asfarray([-100,-100]))
p4=myLattice.addPointMass(1,np.asfarray([-100,100]))
p5=myLattice.addPointMass(1,np.asfarray([0,0]))
myLattice.connect(p1,p2,4)
myLattice.connect(p1,p3,4)
myLattice.connect(p2,p4,4)
myLattice.connect(p3,p4,4)
myLattice.connect(p5,p1,4)
myLattice.connect(p5,p2,4)
myLattice.connect(p5,p3,4)
myLattice.connect(p5,p4,4)




#myLattice.startNewFrame()
#myLattice.addImpulse(p2,np.asfarray([10,0]))
#myLattice.resolveFrame(1)

myLattice.angularVelocity=0.5

myRenderer.addLattice(myLattice)

myRenderer.loop()

t=time.time()
timeMultiplier=4
while True:
  myLattice.startNewFrame()
  #myLattice.addImpulse(p1,np.asfarray([-1,0.1]))
  #myLattice.addImpulse(p2,np.asfarray([-1,-1]))
  #myLattice.addImpulse(p3,np.asfarray([1,1]))
  #myLattice.addImpulse(p4,np.asfarray([1,-1]))
  myLattice.resolveFrame((time.time()-t)*timeMultiplier)
  t=time.time()
  myRenderer.loop()