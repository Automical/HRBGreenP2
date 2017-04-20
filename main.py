#Main code file for servo control through inverse kinematics
#Inverse works (ignoring limitations caused by mechanical issues)
#Mapping coordinates onto the paper does not


#-----------------------------------------------------------------------------------------------
#Imports
from joy import *
import numpy as np
from numpy import asarray
from scipy.linalg import expm as expM
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from time import sleep
import ast
from math import *
from util import rigidFix
#-----------------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------------
#Global Constants
num_motors = 3
num_motors_arm = 3
dist_thresh = 1
dt_max = .25
step = .25
filename = 'h.txt'
#-----------------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------------
#Given Kinematics Code



#-----------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------
#Arm Representation CLass
class Arm( object ):

  def __init__(self):
    # link lengths in cm
    self.l1 = 28
    self.l2 = 35
    self.l3 = 20
    self.l4 = 15
    self.l5 = 47
    self.l6 = 10
    self.l7 = 10
    self.l8 = 5
    
    
  def getToolAng(self, theta1,theta2,theta5):
    #print("theta1,2,5",theta1,theta2,theta5)

    rc = self.l2*cos(theta2) + self.l6
    zc = self.l2*sin(theta2)
    #print("Rc,Zc",rc,zc)
    
    ra = self.l1*cos(theta1) - self.l7
    za = self.l1*sin(theta1)
    #print("Ra,Za",ra,za)
    
    xc = rc
    yc = zc
    
    xa = ra
    ya = za
    #print("Xa,Ya",xa,ya)
    
    d = dist2(xc,yc,xa,ya)
    #print("D:",d)

    a = (self.l3**2-self.l4**2+d**2)/(2*d)
    #print("A:",a)
    
    h = sqrt(self.l3**2-a**2)
    
    x2 = xa +a*(xc-xa)/d
    y2 = ya +a*(yc-ya)/d
    
    x31 = x2 + h * (yc-ya)/d
    y31 = y2 - h * (xc-xa)/d
    
    x32 = x2 - h * (yc-ya)/d
    y32 = y2 + h * (xc-xa)/d
    
    if (y31 > y32):
      yb = y31
      xb = x31
    else:
      yb = y32
      xb = x32
      
    theta4 = atan2(yc-yb,xc-xb)
    
    r = xb + self.l5*cos(theta4)
    z = yb + self.l5*sin(theta4)
    
    x = r*cos(theta5)
    y = r*sin(theta5)
    
    tool = asarray([[x],[y],[z]])
    return tool
    
  def getTool(self,ang):
    rc = self.l2*cos(ang[1]) + self.l6
    zc = self.l2*sin(ang[1])
    re = self.l5*cos(ang[3])+rc
    ze = self.l5*sin(ang[3])+zc
  
    
  
    toolX = re*cos(ang[4])
    toolY = re*sin(ang[4])
    
    toolZ = ze
    
    return asarray([[toolX],[toolY],[toolZ]])
  
  def angFromEnd(self,x,y,z):
    ang = [90,90,90,90,90];
    
    r = dist2(0,0,x,y)
    
    theta5 = atan2(y,x)
    ang[4]=theta5
    x=r
    
    d1 = dist2(self.l6,0,x,z)
    #print("d1 = %f" % d1)
    beta = acos((self.l2**2+d1**2-self.l5**2)/(2*self.l2*d1))
    #print("beta = %f" % degrees(beta))
    alpha = atan2(z,x-self.l6)
    #print("alpha = %f" % degrees(alpha))
    theta2 = alpha + beta;
    #print("theta2 = %f" % degrees(theta2))

    
    ang[1] = theta2
    
    xc = self.l2*cos(theta2) + self.l6
    zc = self.l2*sin(theta2)
    
    #print("xc = %f" % xc)
    #print("zc = %f" % zc)
    
    theta4 = -atan2(zc-z,x-xc)
    ang[3]=theta4;
    

    #print("theta4 = %f" % degrees(theta4))
    xb = xc - self.l4*cos(theta4)
    zb = zc - self.l4*sin(theta4)

    #print("xb = %f" % xb)
    #print("zb = %f" % zb)
    d2 = dist2(-self.l7,0,xb,zb)
    #print("d2 = %f" % d2)
    
    
    
    gamma = acos((self.l1**2+d2**2-self.l3**2)/(2*self.l1*d2))
    #print("gamma = %f" % degrees(gamma))
    delta = atan2(zb,xb-(-self.l7))
    #print("delta = %f" % degrees(delta))
    theta1 = gamma+delta
    
    if (theta1<(pi/2)):
      theta1 = (pi/2-theta1) + pi/2
    
    #print("theta1= %f" % degrees(theta1))
    ang[0] = theta1
    
    xa = self.l1*cos(theta1) - self.l7
    za = self.l1*sin(theta1)
    
    theta3 = atan2(zb-za,xb-xa)
    #print("theta3= %f" % degrees(theta3))

    ang[2] = theta3
    
    
    return ang
  
  def plot3D(self,ang):
    global ax
    rc = self.l2*cos(ang[1]) + self.l6
    zc = self.l2*sin(ang[1])
    re = self.l5*cos(ang[3])+rc
    ze = self.l5*sin(ang[3])+zc
    ra = self.l1*cos(ang[0]) - self.l7
    za = self.l1*sin(ang[0])
    rb = self.l3*cos(ang[2])+ra
    zb = self.l3*sin(ang[2])+za
    
    xc = rc*cos(ang[4])
    yc = rc*sin(ang[4])
    
    xe = re*cos(ang[4])
    ye = re*sin(ang[4])
    
    xa = ra*cos(ang[4])
    ya = ra*sin(ang[4])
    
    xb = rb*cos(ang[4])
    yb = rb*sin(ang[4])
    
    
    x = [-self.l7*cos(ang[4]), xa]
    z = [0,za]
    y = [-self.l7*sin(ang[4]), ya]
    ax.plot(x,y,z,label='l1')
    x = [xa, xb]
    y = [ya,yb]
    z = [za, zb]
    ax.plot(x,y,z,label='l3')
    x = [xb,xc]
    y = [yb,yc]
    z = [zb,zc]
    ax.plot(x,y,z,label='l4')
    x = [self.l6*cos(ang[4]), xc]
    y = [self.l6*sin(ang[4]), yc]
    z = [0,zc]
    ax.plot(x,y,z,label='l2')
    x = [xc,xe]
    y = [yc,ye]
    z = [zc,ze]
    ax.plot(x,y,z,label='l5')
  
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
    translate = t[0:3,3:4]
    
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
  
  def convertPoint2(self,x,y,z):
    v1 = self.p3 - self.p1
    v2 = self.p2 - self.p1
    
    cp = np.cross(v1.T, v2.T)
    #print(cp)
    #a, b, c = cp
    a = cp[0][0]
    b = cp[0][1]
    c = cp[0][2]
    d = np.dot(cp, self.p3)
    #Z = (d - a * X - b * Y) / c
    #d = cz+ax+by
    vector_norm = a*a + b*b + c*c
    
    normal_vector = np.array([a, b, c]) / np.sqrt(vector_norm)
    point_in_plane = np.array([a, b, c]) / vector_norm

    points = np.column_stack((x, y, z))
    points_from_point_in_plane = points - point_in_plane
    proj_onto_normal_vector = np.dot(points_from_point_in_plane,
                                     normal_vector)
    proj_onto_plane = (points_from_point_in_plane -
                       proj_onto_normal_vector[:, None]*normal_vector)

    return point_in_plane + proj_onto_plane
    #x=x-self.p1[0]
    #y=y-self.p1[1]
    #z=z-self.p1[2]
    v = asarray([[x],[y],[z]])
    #x_out = self.basis1*v
    #y_out = self.basis2*v
    #return [x,y]

#-----------------------------------------------------------------------------------------------
    
    
class Controller(object):
  def __init__(self,app):
    self.app = app
  
  #return change in angle for given arm position, angle, and endpoint
  #todo: dynamic dt selection
  '''
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
    '''
  def generateToolDelta(self,arm,ang,end):
    tool = arm.getTool(ang)
    diff = end-tool
    if (abs(diff[0])>.2):
      if (diff[0])>0:
        tx = .2
      else:
        tx=-.2
    else:
      tx = 0
    if (abs(diff[1])>.2):
      if (diff[1])>0:
        ty = .2
      else:
        ty=-.2
    else:
      ty = 0
      
    if (abs(diff[2])>.2):
      if (diff[2])>0:
        tz = .2
      else:
        tz=-.2
    else:
      tz = 0
  
    return asarray([[tx],[ty],[tz]])
  
def dist(tool,end):
  return sqrt((tool[0]-end[0])**2+(tool[1]-end[1])**2+(tool[2]-end[2])**2)

def dist2(x1,y1,x2,y2):
  return sqrt((x2-x1)**2+(y2-y1)**2);

#-----------------------------------------------------------------------------------------------
  
  
#-----------------------------------------------------------------------------------------------
#Plan to drive robot arm to point in relative 3D space
class GoToPoint(Plan):
  def __init__(self,app,*arg,**kw):
    Plan.__init__(self,app,*arg,**kw)
    self.app = app
    self.arm = self.app.arm
    self.ang = self.app.ang
    self.control = Controller(self.app)
    self.setEnd(self.arm.getTool(self.ang))
  
  def setEnd(self,end):
    self.end = end
    print("Set target endpoint to:",self.end)
    
  def behavior(self):
    dist_to = dist(self.arm.getTool(self.ang),self.end)
    
    while(dist_to > dist_thresh):
      dist_to = dist(self.arm.getTool(self.ang),self.end)
      #da = self.control.generateAngleDelta(self.arm,self.ang,self.end)
      #self.ang = self.ang + da
      dpos = self.control.generateToolDelta(self.arm,self.ang,self.end)
      tool  = self.arm.getTool(self.ang)
      tool2 = tool+dpos
      self.ang = self.arm.angFromEnd(tool2[0],tool2[1],tool2[2])
      
      print("Setting0:",self.ang[0]," Setting1:",self.ang[1]," Setting4:",self.ang[4])
      self.app.ser[0].set_pos(worldAngtoRobAng(degrees(self.ang[0])))
      self.app.ser[1].set_pos(worldAngtoRobAng(degrees(self.ang[1])))
      self.app.ser[2].set_pos(-worldAngtoRobAng(degrees(self.ang[4])))
      
      #for i in range(0,num_motors_arm):
      #  self.app.ser[i].set_pos(self.ang[i])
        
      yield self.forDuration(.001)

#-----------------------------------------------------------------------------------------------
class DrawStrokes(Plan):
  def __init__(self,app,*arg,**kw):
    Plan.__init__(self,app,*arg,**kw)
    self.app = app
    self.point_plan = GoToPoint(self.app)
    self.strokes = [[]]

  def setStrokes(self,s):
    self.strokes = s

  def behavior(self):
    print(self.strokes)
    print(self.point_plan)
    
    for s in self.strokes:
      for p in s:
        print("Going to",p)
        self.point_plan.setEnd(self.app.paper.convertPoint(p[0],p[1]))
        yield self.point_plan

      
#-----------------------------------------------------------------------------------------------
#Rigid Body Definition

#Numerical Solution to Wahba's problem
#p1 is origin of measured paper, others are around paper clockwise
#r1-r4 are original, corresponding paper points
#returns 4x4 rigid body transform
'''
def calcRigidTransform(p1,p2,p3,p4,r1,r2,r3,r4):
  B = np.matmul(p1,r1.T)+np.matmul(p2,r2.T)+np.matmul(p3,r3.T)+np.matmul(p4,r4.T)
  U, s, V = np.linalg.svd(B, full_matrices=True)
  M = np.matrix([[1,0,0],[0,1,0],[0,0,np.linalg.det(U)*np.linalg.det(V)]])
  R = np.matmul(U,np.matmul(M,V))
  
  Ro = np.zeros((4,4))
  Ro[0:3,0:3] = R
  Ro[0:4,3:4] = np.matrix([[p1[0]],[p1[1]],[p1[2]],[1]])
  
  #print("Ro",Ro,"p1",p1,"p1r",np.matmul(Ro,p1),"p2",p2,"p2r",np.matmul(Ro,p2),"p3",p3,"p3r",np.matmul(Ro,p3),"p4",p4,"p4r",np.matmul(Ro,p4))
  
  print(Ro)
  
  return Ro
'''

def calcRigidTransform(p1,p2,p3,p4,r1,r2,r3,r4):
  P = [p1,p2,p3,p4]
  R = [r1,r2,r3,r4]


def getToolLoc(arm,ser):
  ang = np.zeros(num_motors_arm)
  for i in range(0,num_motors_arm):
    ang[i] = ser[i].get_pos()
   

  for i in range(0,len(ang)):
    ang[i] = radians(robAngToWorldAng(ang[i]))
    
  
  ang = [ang[0], ang[1], 0, 0, ang[2]]
  
  
  
  
  return arm.getToolAng(ang[0],ang[1],ang[4])
       
   
   
#Robot servo angles may not 1-1 correspond to arm frame angles
def robAngToWorldAng(ang):
  if ang == 0:
    return 90
  return -(ang/100.0) + 90
  
def worldAngtoRobAng(ang):
  
  return int((ang - 90)*-100)


#plots workspace
def plot_cuboid(center, size):
  global ax
  """
     Create a data array for cuboid plotting.


     ============= ================================================
     Argument      Description
     ============= ================================================
     center        center of the cuboid, triple
     size          size of the cuboid, triple, (x_length,y_width,z_height)
     :type size: tuple, numpy.array, list
     :param size: size of the cuboid, triple, (x_length,y_width,z_height)
     :type center: tuple, numpy.array, list
     :param center: center of the cuboid, triple, (x,y,z)
 """
  # suppose axis direction: x: to left; y: to inside; z: to upper
  # get the (left, outside, bottom) point
  ox, oy, oz = center
  l, w, h = size

  x = np.linspace(ox-l/2,ox+l/2,num=2)
  y = np.linspace(oy-w/2,oy+w/2,num=2)
  z = np.linspace(oz-h/2,oz+h/2,num=2)
  x1, z1 = np.meshgrid(x, z)
  y11 = np.ones_like(x1)*(oy-w/2)
  y12 = np.ones_like(x1)*(oy+w/2)
  x2, y2 = np.meshgrid(x, y)
  z21 = np.ones_like(x2)*(oz-h/2)
  z22 = np.ones_like(x2)*(oz+h/2)
  y3, z3 = np.meshgrid(y, z)
  x31 = np.ones_like(y3)*(ox-l/2)
  x32 = np.ones_like(y3)*(ox+l/2)

  # outside surface
  ax.plot_wireframe(x1, y11, z1, color='y', rstride=1, cstride=1, alpha=0.6)
  # inside surface
  ax.plot_wireframe(x1, y12, z1, color='y', rstride=1, cstride=1, alpha=0.6)
  # bottom surface
  ax.plot_wireframe(x2, y2, z21, color='y', rstride=1, cstride=1, alpha=0.6)
  # upper surface
  ax.plot_wireframe(x2, y2, z22, color='y', rstride=1, cstride=1, alpha=0.6)
  # left surface
  ax.plot_wireframe(x31, y3, z3, color='y', rstride=1, cstride=1, alpha=0.6)
  # right surface
  ax.plot_wireframe(x32, y3, z3, color='y', rstride=1, cstride=1, alpha=0.6)



class DrawPlan( Plan ):
  def __init__(self,app,*arg,**kw):
    global ax
    Plan.__init__(self,app,*arg,**kw)
    fig = gcf()
    ax = fig.gca(projection='3d')
    plt.show()
  def behavior(self):
    global ax
    # Clear the screen
    ax.clear()

    # Draw workspace
    plot_cuboid((19.24,0,15.24),(30.48,30.48,30.48))
    
    self.app.arm.plot3D(self.app.ang)
    
    x = asarray([float(self.app.paper.p1[0]),float(self.app.paper.p2[0]),float(self.app.paper.p3[0]),float(self.app.paper.p4[0]),float(self.app.paper.p1[0])])
    y = asarray([float(self.app.paper.p1[1]),float(self.app.paper.p2[1]),float(self.app.paper.p3[1]),float(self.app.paper.p4[1]),float(self.app.paper.p1[1])])
    z = asarray([float(self.app.paper.p1[2]),float(self.app.paper.p2[2]),float(self.app.paper.p3[2]),float(self.app.paper.p4[2]),float(self.app.paper.p1[2])])
    
    plt.draw()
    plt.pause(.0001)
    show()
    yield(self.forDuration(.1))
    

#-----------------------------------------------------------------------------------------------
#App
class GreenApp( JoyApp ):
  def __init__(self,*arg,**kw):
    JoyApp.__init__( self, confPath="$/cfg/JoyApp.yml", *arg, **kw) 

    #self.ang = [0,pi/4,pi/4]
    self.ser_t = self.robot.at
    #Replace with real motor value
    
    self.arm = Arm()
    xi = 30
    yi = 0
    zi = 35
  
    
    self.ang = self.arm.angFromEnd(xi,yi,zi)
    self.paper = Paper()
    #Order in the direction of joints
    for x in self.ser_t:
      print(x)
    self.ser = [self.ser_t.K21,self.ser_t.Nx17,self.ser_t.Nx21];
    self.ang[0] = radians(robAngToWorldAng(self.ser_t.K21.get_pos()))
    self.ang[1] = radians(robAngToWorldAng(self.ser_t.Nx17.get_pos()))
    self.ang[4] = radians(robAngToWorldAng(self.ser_t.Nx21.get_pos()))
    #self.ser = [0,0,0];
    
    self.p1_d = np.matrix([[0],[0],[0]])
    self.p2_d = np.matrix([[20.32],[0],[0]])
    self.p3_d = np.matrix([[20.32],[27.94],[0]])
    self.p4_d = np.matrix([[0],[27.94],[0]])
    self.p1 = self.p1_d
    self.p2 = self.p2_d
    self.p3 = self.p3_d
    self.p4 = self.p4_d
    
    self.T = np.matrix([[ 3./5, 4./5, 0., 50.],[-4./5, 3./5, 0., 0.],[0.,0.,1.,0.],[0.,0.,0.,1.]])
    self.points_on = [[],[]]
    self.points_off = [[],[]]
    
    self.point_plan = GoToPoint(self)
    self.stroke_plan = DrawStrokes(self)
    self.draw_plan = DrawPlan(self)
    
  def onStart( self ):
    print("start")
    self.draw_plan.start()
    
  def onEvent( self, evt ):
    if evt.type == KEYDOWN:
      if evt.key == K_UP:
        print("up")
      elif evt.key == K_o:
        for i in range(0,num_motors):
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
        self.T = calcRigidTransform(self.p1,self.p2,self.p3,self.p4,self.p1_d,self.p2_d,self.p3_d,self.p4_d)
        self.paper.update_paper(self.T)
      elif evt.key == K_f:
        progress("read file")
        f = open(filename, 'r')
        l = f.readlines()[0]
        l=l.replace('(','[')
        l=l.replace(')',']')
        l=l.replace(' ','')
        z=ast.literal_eval(l)
        self.points_on = z;
        self.stroke_plan.setStrokes(self.points_on)
      elif evt.key == K_k:
        if self.point_plan.isRunning() == False:
          d = input("target as list")
          print("Going to some rando point!")
          self.point_plan.setEnd(asarray(d))
          self.point_plan.start()
        else:
          print("Point plan already running")
      elif evt.key == K_l:
        #pass
        print("Trying to draw points")
        #self.point_plan.stop()
        self.stroke_plan.start()
      elif evt.key == K_ESCAPE:
        self.stop()
      elif evt.key==K_b:
        print(getToolLoc(self.arm,self.ser))
      elif evt.key==K_v:
        #z = getToolLoc(self.arm,self.ser)
        #print(self.paper.convertPoint2(z[0],z[1],z[2]))
        for s in self.points_on:
              for p in s:
                print("printing point")
                print("Going to ",p)
                print("Normal Coords ", self.paper.convertPoint(p[0],p[1]))
                





if __name__=="__main__":
  
  import sys
  app=GreenApp(robot = dict(count=num_motors
                   ,port=dict(TYPE='TTY', glob="/dev/ttyACM*", baudrate=115200)))
  app.run()