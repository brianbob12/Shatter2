
import time
import numpy as np
import Renderer
import Engine

myRenderer=Renderer.Renderer(1024,1024)

myLattice,edgePoints=Engine.rect(0,0,200,400,10,5e-3,40,1)

myRenderer.addLattice(myLattice)

chosen=0


myRenderer.loop()

myLattice.startNewFrame()
myLattice.addImpulse(edgePoints[chosen],np.asfarray([-100,0]))
myLattice.resolveFrame(1)

t=time.time()
timeMultiplier=0.2
while True:
  myLattice.startNewFrame()
  myLattice.addImpulse(edgePoints[chosen],np.asfarray([-10,0]))
  myLattice.resolveFrame((time.time()-t)*timeMultiplier)
  t=time.time()
  myRenderer.loop()
  myRenderer.drawPoint(edgePoints[chosen],myLattice,(255,0,255))
  myRenderer.update()
  