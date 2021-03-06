
import time
import numpy as np
import Renderer
import Engine

myRenderer=Renderer.Renderer(1024,1024)
mySystem=Engine.System(timeMultiplier=0.1,g=0)


myLattice,edgePoints=Engine.regularRect(300,0,200,800,10,1e-2,1)
myRenderer.addLattice(myLattice)
mySystem.addLattice(myLattice)


myCircle,circlePerimeter=Engine.circle(-200,-300,50,10,2e-2,20,1)

myRenderer.addLattice(myCircle)
mySystem.addLattice(myCircle)

myCircle.velocity=np.asfarray([20,0])
myLattice.angle=3.141592/2

t=time.time()
i=0
while True:
  tp=time.time()
  #mySystem.update(tp-t)
  mySystem.update(1)
  t=tp
  myRenderer.loop()
  #myRenderer.saveFrame(f"tests/B1/{i:04d}.png")
  """
  for j in range(min(i,len(myLattice.getSurfaceConnections()))):
    myRenderer.drawPoint(myLattice.getSurfaceConnections()[j].p1,myLattice,(255,0,0))
    myRenderer.drawPoint(myLattice.getSurfaceConnections()[j].p2,myLattice,(255,0,0))
  myRenderer.update()
  """
  i+=1
  