
import time
import numpy as np
import Renderer
import Engine

myRenderer=Renderer.Renderer(1024,1024)
mySystem=Engine.System(timeMultiplier=0.05,g=0)


myLattice,edgePoints=Engine.rect(300,0,200,400,10,1e-2,25,1)
myRenderer.addLattice(myLattice)
mySystem.addLattice(myLattice)


myCircle,circlePerimeter=Engine.circle(-200,0,50,10,2e-2,20,1)

myRenderer.addLattice(myCircle)
mySystem.addLattice(myCircle)

myCircle.velocity=np.asfarray([10,0])
myLattice.angle=3.141592/4

t=time.time()
i=0
while True:
  tp=time.time()
  #mySystem.update(tp-t)
  mySystem.update(1)
  t=tp
  myRenderer.loop()
  myRenderer.saveFrame(f"tests/B0/{i:04d}.png")
  """
  for j in range(min(i,len(myLattice.getSurfaceConnections()))):
    myRenderer.drawPoint(myLattice.getSurfaceConnections()[j].p1,myLattice,(255,0,0))
    myRenderer.drawPoint(myLattice.getSurfaceConnections()[j].p2,myLattice,(255,0,0))
  myRenderer.update()
  """
  i+=1
  