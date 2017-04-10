#-----------------------------------------------------------------------------------------------
#Imports
from joy import *
import numpy as np
from scipy.linalg import expm as expM
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from time import sleep
import ast
#-----------------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------------
#Global Constants
num_motors = 3
num_arm_motors = 3
dist_thresh = 1
dt_max = .25

#-----------------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------------
#Given Kinematics Code

def seToSE( x ):
  """
  Convert a twist (a rigid velocity, element of se(3)) to a rigid
  motion (an element of SE(3))
  
  INPUT:
    x -- 6 sequence
  OUTPUT:
    result -- 4 x 4  

  """
  x = asarray(x,dtype=float)
  if x.shape != (6,):
    raise ValueError("shape must be (6,); got %s" % str(x.shape))
  #
  return expM(screw(x))

def screw( v ):
  """
  Convert a 6-vector to a screw matrix 
  
  The function is vectorized, such that:
  INPUT:
    v -- N... x 6 -- input vectors
  OUTPUT:
    N... x 4 x 4  
  """
  v = asarray(v)
  z = zeros_like(v[0,...])
  return array([
      [ z, -v[...,5], v[...,4], v[...,0] ],
      [ v[...,5],  z,-v[...,3], v[...,1] ],
      [-v[...,4],  v[...,3], z, v[...,2] ],
      [ z,         z,        z, z] ])

def unscrew( S ):
  """
  Convert a screw matrix to a 6-vector
  
  The function is vectorized, such that:
  INPUT:
    S -- N... x 4 x 4 -- input screws
  OUTPUT:
    N... x 6
  
  This is the "safe" function -- it tests for screwness first.
  Use unscrew_UNSAFE(S) to skip this check
  """
  S = asarray(S)
  assert allclose(S[...,:3,:3].transpose(0,1),-S[...,:3,:3]),"S[...,:3,:3] is skew"
  assert allclose(S[...,3,:],0),"Bottom row is 0"
  return unscrew_UNSAFE(S)

def jacobian_cdas( func, scl, lint=0.8, tol=1e-12, eps = 1e-30, withScl = False ):
  """Compute Jacobian of a function based on auto-scaled central differences.
  
  INPUTS:
    func -- callable -- K-vector valued function of a D-dimensional vector
    scl -- D -- vector of maximal scales allowed for central differences
    lint -- float -- linearity threshold, in range 0 to 1. 0 disables
         auto-scaling; 1 requires completely linear behavior from func
    tol -- float -- minimal step allowed
    eps -- float -- infinitesimal; must be much smaller than smallest change in
         func over a change of tol in the domain.
    withScl -- bool -- return scales together with Jacobian
  
  OUTPUTS: jacobian function 
    jFun: x --> J (for withScale=False)
    jFun: x --> J,s (for withScale=True)
    
    x -- D -- input point
    J -- K x D -- Jacobian of func at x
    s -- D -- scales at which Jacobian holds around x
  """
  scl = abs(asarray(scl).flatten())
  N = len(scl)  
  lint = abs(lint)
  def centDiffJacAutoScl( arg ):
    """
    Algorithm: use the value of the function at the center point
      to test linearity of the function. Linearity is tested by 
      taking dy+ and dy- for each dx, and ensuring that they
      satisfy lint<|dy+|/|dy-|<1/lint
    """
    x0 = asarray(arg).flatten()    
    y0 = func(x0)
    s = scl.copy()
    #print "Jac at ",x0
    idx = slice(None)
    dyp = empty((len(s),len(y0)),x0.dtype)
    dyn = empty_like(dyp)
    while True:
      #print "Jac iter ",s
      d0 = diag(s)
      dyp[idx,:] = [ func(x0+dx)-y0 for dx in d0[idx,:] ]
      dypc = dyp.conj()
      dyn[idx,:] = [ func(x0-dx)-y0 for dx in d0[idx,:] ]
      dync = dyn.conj()      
      dp = sum(dyp * dypc,axis=1)
      dn = sum(dyn * dync,axis=1)
      nul = (dp == 0) | (dn == 0)
      if any(nul):
        s[nul] *= 1.5
        continue
      rat = dp/(dn+eps)
      nl = ((rat<lint) | (rat>(1.0/lint)))
      # If no linearity violations found --> done
      if ~any(nl):
        break
      # otherwise -- decrease steps
      idx, = nl.flatten().nonzero()
      s[idx] *= 0.75
      # Don't allow steps smaller than tol
      s[idx[s[idx]<tol]] = tol
      if all(s[idx]<tol):
        break
    res = ((dyp-dyn)/(2*s[:,newaxis])).T
    if withScl:
      return res, s
    return res
  return centDiffJacAutoScl 

#-----------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------
#Arm Representation CLass
class Arm( object ):
  """
  class Arm
  
  Represents a series manipulator made of several segments.
  Each segment is graphically represented by a wireframe model
  """
  def __init__(self):
    # link lengths
    self.ll = asarray([20,20,20])
    d = 0.2

    self.geom = [( asarray([[0,0,0,1]]) ).T ]
    #
    # Build twist matrices 
    # Build wireframe of all segments in the world coordinates
    #
    tw = []
    LL = 0
    
    orientations = [[0,0,1],[0,1,0],[0,1,0]]
    counter = 0
    for n,ll in enumerate(self.ll):
      # Scale the geometry to the specifies link length (ll)
      # Shift it to the correct location (LL, sum of the ll values)
      self.geom.append([( asarray([[0,0,LL,1]]) ).T ])

      # Compute the twist for this segment; 
      # twists alternate between the z and y axes
      w = orientations[counter]
      # Velocity induced at the origin
      v = -cross(w,[0,0,LL])
      # Collect the twists
      tw.append( concatenate([v,w],0) )
      # Accumulate the distance along the arm
      LL += ll

      counter = counter + 1
    # Build an array of collected twists
    self.tw = asarray(tw)
    self.tool = asarray([0,0,LL,1]).T
    
    
    # overwrite method with jacobian function
    self.getToolJac = jacobian_cdas( 
      self.getTool, ones(self.tw.shape[0])*0.05 
    )
  
  def at( self, ang ):
    """
    Compute the rigid transformations for a multi-segment arm
    at the specified angles
    """
    ang = asarray(ang)[:,newaxis]
    tw = ang * self.tw
    
    A = [identity(4)]
    for twi in tw:
      M = seToSE(twi)
      A.append(dot(A[-1],M))
    return A
    
  def getTool( self, ang ):
    """
    Get "tool tip" position in world coordinates
    """
    # Get the rigid transformation for the last segment of the arm
    M = self.at(ang)[-1]
    return dot(M, self.tool)
  
  def getToolJac( self, ang ):
    """
    Get "tool tip" Jacobian by numerical approximation
    
    NOTE: implementation is a placeholder. This method is overwritten
    dynamically by __init__() to point to a jacobian_cdas() function
    """
    raise RuntimeError("uninitialized method called")
    
  def plot3D( self, ang ):
    global ax
    """
    Plot arm in 3 views
    """
    x = []
    y = []
    z = []
    A = self.at(ang)
    for a,g in zip(A, self.geom):
      ng = dot(a,g)
      x.append(float(ng[0]))
      y.append(float(ng[1]))
      z.append(float(ng[2]))
    tp = self.getTool(ang)
    x.append(float(tp[0]))
    y.append(float(tp[1]))
    z.append(float(tp[2]))
    ax.plot( x, y, z,  zdir='z' )
    ax.plot( x, y, z, 'hk', zdir='z' )
    ax.set_xlim(-30,30)
    ax.set_xlabel('x')
    ax.set_ylim(-30,30)
    ax.set_ylabel('y')
    ax.set_zlim(-30,30)
    ax.set_zlabel('z')
    tx = asarray([0,(tp[0])])
    ty = asarray([0,(tp[1])])
    tz = asarray([0,(tp[2])])
    ax.plot(tx,ty,tz,'hr',zdir='z')
    
#-----------------------------------------------------------------------------------------------
#Representation of paper in 3D space
#Conversions from 2D paper space to 3D arm space
class Paper(object):
  def __init__(self):
    
    #Generate default paper coordinates
    self.p1_d = np.matrix([[0],[0],[0]])
    self.p2_d = np.matrix([[20.32],[0],[0]])
    self.p3_d = np.matrix([[20.32],[27.94],[0]])
    self.p4_d = np.matrix([[0],[27.94],[0]])
    
    self.p1 = self.p1_d
    self.p2 = self.p2_d
    self.p3 = self.p3_d
    self.p4 = self.p4_d
    
    self.update_basis()
  
  #Converts representation of paper to rotated and translated paper in 3D space
  #t - rotation/translation in homogenous coordinates
  def update_paper(self,t):
    rotate = t[0:3,0:3]
    translate = t[0:3,3:4]/10.0
    
    self.p1 = rotate*self.p1_d + translate
    self.p2 = rotate*self.p2_d + translate
    self.p3 = rotate*self.p3_d + translate
    self.p4 = rotate*self.p4_d + translate
    
    self.update_basis()
    
  def update_basis(self):
    self.basis1 = self.p2-self.p1
    self.basis1 = self.basis1/np.sqrt(np.multiply(self.basis1,self.basis1)).sum()
    self.basis2 = self.p4-self.p1
    self.basis2 = self.basis2/np.sqrt(np.multiply(self.basis2,self.basis2)).sum()
  
  #converts a point in 2D space on the paper to a point in 3D space relative to the arm
  def convertPoint(self,x,y):
    return self.p1 + self.basis1*x+self.basis2*y
  
#-----------------------------------------------------------------------------------------------
    
    
class Controller(object):
  def __init__(self,app):
    self.app = app
  
  #return change in angle for given arm position, angle, and endpoint
  #todo: dynamic dt selection
  def generateAngleDelta(self,arm,ang,end):
    tool = arm.getTool(ang)
    Jt = a.getToolJac(ang)
    
    dt = .10
    
    dx = 0
    dy = 0
    dz = 0
    
    if (end[0]-tool[0])>0:
      dx = dt
    else:
      dx = -dt

    if (end[1]-tool[1])>0:
      dy = dt
    else:
      dy = -dt
    
    if (end[2]-tool[2])>0:
      dz = dt
    else:
      dz = -dt
      
    d = asarray([tx,ty,tz])
    
    da = dot(pinv(Jt)[:,:len(d)],d)
    
    return da
    
  def dist(tool,end):
    return sqrt((tool[0]-end[0])**2+(tool[1]-end[1])**2+(tool[2]-end[2])**2)

  
#-----------------------------------------------------------------------------------------------
  
  
#-----------------------------------------------------------------------------------------------
#Plan to drive robot arm to point in relative 3D space
class GotoPoint(Plan):
  def __init__(self,app,*arg,**kw):
    Plan.__init__(self,app,*arg,**kw)
    self.app = app
    self.arm = self.app.arm
    self.ang = self.app.ang
    self.control = Controller(self.app)
    self.setEnd(self.arm.getTool(self.ang))
  
  def setEnd(self,end):
    self.end = end
    
  def behavior(self):
    dist = self.control.dist(self.arm.getTool(self.ang),self.end)
    
    while(dist > dist_thresh):
      dist = self.control.dist(self.arm.getTool(self.ang),self.end)
      da = self.control.generateAngleDelta(self.arm,self.ang,self.end)
      self.ang = self.ang + da
      
      for (i in range(1,num_motors_arm)):
        self.app.ser[i].set_pos(self.ang[i])
        
      yield self.forDuration(.1)

      
#-----------------------------------------------------------------------------------------------
#Rigid Body Definition

#Numerical Solution to Wahba's problem
#p1 is origin of measured paper, others are around paper clockwise
#r1-r4 are original, corresponding paper points
#returns 4x4 rigid body transform
def calcRigidTransform(p1,p2,p3,p4,r1,r2,r3,r4):
  B = np.matmul(p1,r1.T)+np.matmul(p2,r2.T)+np.matmul(p3,r3.T)+np.matmul(p4,r4.T)
  U, s, V = np.linalg.svd(B, full_matrices=True)
  M = np.matrix([[1,0,0],[0,1,0],[0,0,np.linalg.det(U)*np.linalg.det(V)]])
  R = np.matmul(U,np.matmul(M,V))
  
  Ro = np.zeros(4,4)
  Ro[0:3,0:3] = R
  Ro[0:3,3:4] = np.matrix([[p1[0]],[p1[1]],[p1[2]],[1]])
  return Ro
  
def getToolLoc(arm,ser):
  ang = np.zeros(1,num_motors_arm)
  for (i in range(1,num_motors_arm):
     ang[i] = ser[i].get_pos()
   
   ang = robAngToWorldAng(ang)
   
   returm arm.getTool(ang)
       
   
   
#Robot servo angles may not 1-1 correspond to arm frame angles
def robAngToWorldAng(ang):
     
   return ang/1000.0;
  
#-----------------------------------------------------------------------------------------------
#App
class BlueApp( JoyApp ):
  def __init__(self,*arg,**kw):
    JoyApp.__init__( self, confPath="$/cfg/JoyApp.yml", *arg, **kw) 
    self.ang = [0,pi/4,pi/4]
    self.ser_t = self.robot.at
    #Replace with real motor value
    
    self.arm = Arm()
    
    self.paper = Paper()
    #Order in the direction of joints
    #self.ser = [self.ser_t.Nx45,self.ser_t.Nx23,self.ser_t.Nx4E];
    self.ser = [0,0,0];
    
    self.p1_d = np.matrix([[0],[0],[0]])
    self.p2_d = np.matrix([[20.32],[0],[0]])
    self.p3_d = np.matrix([[20.32],[27.94],[0]])
    self.p4_d = np.matrix([[0],[27.94],[0]])
    self.p1 = self.p1_d
    self.p2 = self.p2_d
    self.p3 = self.p3_d
    self.p4 = self.p4_d
    
    self.T = np.matrix([[ 3./5, 4./5, 0., 50.],[-4./5, 3./5, 0., 0.],[0.,0.,1.,0.],[0.,0.,0.,1.]])

    
  def onStart( self ):
    print("start")
    
  def onEvent( self, evt ):
    print("event")
    
    if evt.type == KEYDOWN:
      if evt.key == K_UP:
        print("up")
      elif evt.key == K_o:
        for (i in range(1,num_motors)):
          self.ser[i].go_slack()
      elif evt.key == K_a:
        progress("Read p1")
        self.p1 = getToolLoc(self.arm,self.ser)
      elif evt.key == K_q:
        progress("Read p2")
        self.p2 = getToolLoc(self.arm,self.ser)
      elif evt.key == K_w:
        progress("Read p3")
        self.p3 = getToolLoc(self.arm,self.ser)
      elif evt.key == K_s:
        progress("Read p4")
        self.p4 = getToolLoc(self.arm,self.ser)
      elif evt.key == K_z:
        self.T = calcRigidTransform(self.p1,self.p2,self.p3,sel.fp4,self.p1_d,self.p2_d,self.p3_d,self.p4_d)
        self.paper.update_paper(self.T)
      elif evt.key == K_f:
        progress("read file")
        








if __name__=="__main__":
  
  import sys
  app=GreenApp(wphAddr=WAYPOINT_HOST, robot = dict(count=num_motors
                   ,port=dict(TYPE='TTY', glob="/dev/ttyACM*", baudrate=115200)))
  app.run()