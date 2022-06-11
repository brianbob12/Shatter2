from turtle import position
from typing import List
from numpy import ndarray,exp
import pygame
import sys
from .Colors import *
from Engine import Lattice, PointMass

class Renderer:
  def __init__(self,width:int,height:int):
    pygame.init()
    self.width:int=width
    self.height:int=height
    self.scale:float=1
    self.xOffset:float=0
    self.yOffset:float=0
    self.screen=pygame.display.set_mode((width,height))
    self.lattices:List[Lattice]=[]

    self.pointRadius=3
    self.connectionThickness=1

    self.dragSpeed:float=1
    self.zoomSpeed:float=0.1

    self.mouseDownLocation:Tuple[int,int]=(0,0)
    self.mouseDown=False

  def clear(self):
    self.screen.fill(black)

  def addLattice(self,lattice:Lattice):
    self.lattices.append(lattice)
  
  def getPixelLocation(self,positionVector:ndarray)-> Tuple[int,int]:
    x=positionVector[0]
    y=positionVector[1]
    x/=self.scale
    y/=self.scale
    x+=self.width/2+self.xOffset
    y+=self.height/2+self.yOffset
    return round(x),round(y)

  def drawLattice(self,lattice:Lattice,color:Tuple[int,int,int]):
    
    for connection in lattice.connections.values():
      end1=self.getPixelLocation(lattice.getAbsolutePosition(connection.p1.relativePos))
      end2=self.getPixelLocation(lattice.getAbsolutePosition(connection.p2.relativePos)) 
      try:
        v=255/(1+exp(-lattice.connectionStresses[connection]))
        compressive=0
        tensile=0
        if lattice.connectionStresses[connection]<0:
          compressive=255-v
        else:
          tensile=v
        c=(255-tensile,255-compressive-tensile,255-compressive)
        try:
          pygame.draw.line(self.screen,c,end1,end2,width=self.connectionThickness)
        except ValueError as e:
          print(e)
          print(v)
      except TypeError as e:
        #numbers are too big
        pass
    
    for pointMass in lattice.points.values():
      self.drawPoint(pointMass,lattice,color)

  def drawPoint(self,point:PointMass,lattice:Lattice,color:Tuple[int,int,int]):
    position=self.getPixelLocation(lattice.getAbsolutePosition(point.relativePos))
    try:
      pygame.draw.circle(self.screen,color,position,self.pointRadius)
    except TypeError as e:
      #numbers are too big
      pass

  def update(self):
    pygame.display.update()
    
  def loop(self):
    for event in pygame.event.get():
      if event.type== pygame.QUIT:
        sys.exit()
      elif event.type == pygame.MOUSEBUTTONDOWN:
        if event.button==3:
          self.mouseDownLocation=pygame.mouse.get_pos()
          self.mouseDown=True
        elif event.button==4:
          #scroll up
          self.scale-=self.zoomSpeed
          #min check
          if self.scale<=0:
            self.scale=self.zoomSpeed
        elif event.button==5:
          #scroll down
          self.scale+=self.zoomSpeed
      elif event.type==pygame.MOUSEBUTTONUP:
        if event.button==3:
          mousePos=pygame.mouse.get_pos()
          self.xOffset-=self.dragSpeed*(self.mouseDownLocation[0]-mousePos[0])
          self.yOffset-=self.dragSpeed*(self.mouseDownLocation[1]-mousePos[1])
          self.mouseDown=False


    if self.mouseDown:
      mousePos=pygame.mouse.get_pos()
      oldXOffset=self.xOffset
      oldYOffset=self.yOffset
      self.xOffset-=self.dragSpeed*(self.mouseDownLocation[0]-mousePos[0])
      self.yOffset-=self.dragSpeed*(self.mouseDownLocation[1]-mousePos[1])
      
    self.clear()
    for lattice in self.lattices:
      self.drawLattice(lattice,green)
    pygame.display.update()

    if self.mouseDown:
      self.xOffset=oldXOffset
      self.yOffset=oldYOffset

  def saveFrame(self,fileName:str):
    pygame.image.save(self.screen,fileName)