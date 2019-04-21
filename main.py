from direct.showbase.ShowBase import ShowBase
from direct.task import Task
from direct.actor.Actor import Actor
from panda3d.core import *

import numpy as np
import math

thetaDegrees = 30
thetaRadians = thetaDegrees * (math.pi / 180.0)

dt = 0.005
time = 0
vel = [0,0,0]

def rotation_matrix(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.
    """
    axis = np.asarray(axis)
    axis = axis / math.sqrt(np.dot(axis, axis))
    a = math.cos(theta / 2.0)
    b, c, d = -axis * math.sin(theta / 2.0)
    aa, bb, cc, dd = a * a, b * b, c * c, d * d
    bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
    return np.array([[aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
                     [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
                     [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]]) 

def setPos(panda3dObj, vec):
	panda3dObj.setPos(vec[0], vec[1], vec[2])

def getInitRefPointPos(theta):
	return [0,-120 * math.sin(theta), 120 * math.cos(theta)]

def getInitPendelPos(theta):
	return [15,-120 * math.sin(theta), 120 * math.cos(theta)]

pos = getInitPendelPos(thetaRadians)

def getInitCamPos(theta):
	return [0,-200 * math.sin(theta), 200 * math.cos(theta)]

def normalize(v):
    norm=np.linalg.norm(v)
    return v/norm

class MyApp(ShowBase):
    def __init__(self):
        ShowBase.__init__(self)

	plightSun = PointLight('Sun')
	plightSun.setColor(VBase4(0.9,0.9,1.0, 1))
	plnpSun = render.attachNewNode(plightSun)
	plnpSun.setPos(0, -1000, 0)
	render.setLight(plnpSun)

	plightMoon = PointLight('Moon')
	plightMoon.setColor(VBase4(1.0,0.9,0.9, 1))
	plnpMoon = render.attachNewNode(plightMoon)
	plnpMoon.setPos(0, 1000, 0)
	render.setLight(plnpMoon)
 
        # Load the environment model.
        self.earth = self.loader.loadModel("untitled.egg")
        # Reparent the model to render.
        self.earth.reparentTo(self.render)
        # Apply scale and position transforms on the model.
	self.earth.setScale(100)
        self.earth.setPos(0, 0, 0)

	self.cam.setPos(0, -300, 0)
	#setPos(self.cam, getInitCamPos(thetaRadians))
        #self.camera.setHpr(0, 0, thetaDegrees)

        self.refPoint = self.loader.loadModel("untitled.egg")
        self.refPoint.reparentTo(self.render)
	self.refPoint.setScale(1)
	self.refPoint.setColor(1,0,0,1)
	setPos(self.refPoint, getInitRefPointPos(thetaRadians))

        self.pendel = self.loader.loadModel("untitled.egg")
        self.pendel.reparentTo(self.render)
	self.pendel.setScale(2)
	self.pendel.setColor(0,1,0,1)
	setPos(self.pendel, getInitPendelPos(thetaRadians))

        # Add the spinCameraTask procedure to the task manager.
        self.taskMgr.add(self.advanceTime, "AdvanceTimeTask")
 
    # Define a procedure to move the camera.
    def advanceTime(self, task):
	global dt
	global time
	global vel
	global pos

	time = time + dt

	print time

        angleDegrees = task.time
        angleRadians = angleDegrees * (math.pi / 180.0)

	#rotate earth
	self.earth.setHpr(angleDegrees, 0, 0)
	
	newRefPoint = np.dot(rotation_matrix([0,0,1],angleRadians), getInitRefPointPos(thetaRadians))

	dif = pos - newRefPoint
	difNorm = np.linalg.norm(dif)
	
	# gravity
	force = - 100 * normalize(newRefPoint)

	# dedx
	force -= 100000 * (1 - 15/difNorm) * dif

	vel += dt * force

	pos += dt * vel

	setPos(self.refPoint, newRefPoint)
	setPos(self.pendel, pos)

        return Task.cont
 
app = MyApp()
app.run()
