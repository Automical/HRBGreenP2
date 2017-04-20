#Inverse kinematics demo that models our actual used robot
#Give program end-point with nested list ie [[x],[y],[z]] to test
#  inverse kinematics

import numpy as np
#4x4 Rigid body transformation for paper representation
PAPER = np.matrix([[ 1/1.414, 0, -1/1.414, 100.],[0., 1,.0, -100.],[1/1.414,0.,1/1.414,0.],[0.,0.,0.,1.]])

PAPER = np.matrix([[  0.21295503,   0.82962318,  -0.51611581,  34.24058021],[ -0.67166225,   0.50793011,   0.53932998, -16.47151874],[  0.70959141,   0.23180248,   0.66539285,  14.14060594],[  0.,           0.,           0. ,          1.        ]])

PAPER = np.matrix([[  5.17839745e-01,  -7.56234943e-01,  -3.99938383e-01,  -4.43822786e-01],[  5.80268561e-01,   6.54016807e-01,  -4.85335361e-01,   4.89009121e+01],[  6.28593983e-01,   1.92542697e-02,   7.77495259e-01,   1.04948072e+01],[  0.00000000e+00,   0.00000000e+00,   0.00000000e+00,   1.00000000e+00]])


#workspace corners
#WORKSPACE = [[0,0,0],[30.48,0,0],[30.48,30.48,0],[0,30.48,0],[0,0,0],[0,0,30.48],[30.48,0,30.48],[30.48,0,0],[0,0,0],[0,0,30.48],[0,30.48,30.48],[0,30.48,0],[0,0,0]

#list of lists of lists/tuples representing pen strokes

PEN_STROKES = [ [(100,0),(100,200)], [[100,100],[200,100]],[[200,00],[200,200]] ]


from scipy.linalg import expm as expM
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from time import sleep
from itertools import product, combinations

from math import *

class ConvertPage(object):
  def __init__(self, o, b1, b2):
    self.origin = o
    self.basis1 = b1
    self.basis2 = b2
        
  def __call__(self, x, y):
    return self.origin + self.basis1 * x + self.basis2 * y


#plots workspace
def plot_cuboid(center, size):
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



def dist2(x1,y1,x2,y2):
  return sqrt((x2-x1)**2+(y2-y1)**2);

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
    print("d1 = %f" % d1)
    beta = acos((self.l2**2+d1**2-self.l5**2)/(2*self.l2*d1))
    print("beta = %f" % degrees(beta))
    alpha = atan2(z,x-self.l6)
    print("alpha = %f" % degrees(alpha))
    theta2 = alpha + beta;
    print("theta2 = %f" % degrees(theta2))

    
    ang[1] = theta2
    
    xc = self.l2*cos(theta2) + self.l6
    zc = self.l2*sin(theta2)
    
    print("xc = %f" % xc)
    print("zc = %f" % zc)
    
    theta4 = -atan2(zc-z,x-xc)
    ang[3]=theta4;
    

    print("theta4 = %f" % degrees(theta4))
    xb = xc - self.l4*cos(theta4)
    zb = zc - self.l4*sin(theta4)
    #xb = self.l6 + self.l2*cos(theta2) + self.l4*cos(theta4)
    #yb = self.l2*sin(theta2) + self.l4*sin(theta4)
    print("xb = %f" % xb)
    print("zb = %f" % zb)
    d2 = dist2(-self.l7,0,xb,zb)
    print("d2 = %f" % d2)
    
    
    
    gamma = acos((self.l1**2+d2**2-self.l3**2)/(2*self.l1*d2))
    print("gamma = %f" % degrees(gamma))
    delta = atan2(zb,xb-(-self.l7))
    print("delta = %f" % degrees(delta))
    theta1 = gamma+delta
    
    if (theta1<(pi/2)):
      theta1 = (pi/2-theta1) + pi/2
    
    print("theta1= %f" % degrees(theta1))
    ang[0] = theta1
    
    xa = self.l1*cos(theta1) - self.l7
    za = self.l1*sin(theta1)
    
    theta3 = atan2(zb-za,xb-xa)
    print("theta3= %f" % degrees(theta3))
    
    #theta3 = acos((self.l1**2+self.l3**2-d2**2)/(2*self.l1*self.l3))
    
    
    
    #if (theta3>(pi/2)):
    #  theta3 = (theta3-pi/2)
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
    
    
    #print("l1 = %f" % dist2(-self.l7,0,xa,ya))
    #print("l2 = %f" % dist2(self.l6,0,xc,yc))
    #print("l3 = %f" % dist2(xa,ya,xb,yb))
    #print("l4 = %f" % dist2(xb,yb,xc,yc))
    #print("l5 = %f" % dist2(xc,yc,xe,ye))
    
def dist(tool,end):
  return sqrt((tool[0]-end[0])**2+(tool[1]-end[1])**2+(tool[2]-end[2])**2)
    
  
def goToPoint(a,ang,end,tip_points, store_points):
  global ax,x,y,z
  print("goto")
  print(ang)
  print(end)
  go = 1
  while (go):

    tool = a.getTool(ang)
    
    diff = end-tool;
    print("goto")
    print("angle")
    print(ang)
    print("end")
    print(end)
    print("diff")
    print(diff)
    if (abs(diff[0])>.1):
      if (diff[0])>0:
        tx = .1
      else:
        tx=-.1
    else:
      tx = 0
    if (abs(diff[1])>.1):
      if (diff[1])>0:
        ty = .1
      else:
        ty=-.1
    else:
      ty = 0
      
    if (abs(diff[2])>.1):
      if (diff[2])>0:
        tz = .1
      else:
        tz=-.1
    else:
      tz = 0
    
    '''
    if (end[1]-tool[1])>0:
      ty = .25
    else:
      ty=-.25
    if (end[2]-tool[2])>0:
      tz = .25
    else:
      tz=-.25
    '''
    
    tool2 = tool + asarray([[tx],[ty],[tz]])
    
    print("tool")
    print(tool)
    print("updated tool")
    print(tool2)
    
    ang = a.angFromEnd(tool2[0],tool2[1],tool2[2])
    '''
    Jt = a.getToolJac(ang)
    if (end[0]-tool[0])>0:
      tx = .25
    else:
      tx=-.25

    if (end[1]-tool[1])>0:
      ty = .25
    else:
      ty=-.25
    if (end[2]-tool[2])>0:
      tz = .25
    else:
      tz=-.25

    d = asarray([tx,ty,tz])

    ang = ang + dot(pinv(Jt)[:,:len(d)],d)
    iteration(a, ang, tip_points, store_points)
    '''
    
    iteration(a,ang,tip_points,store_points)
    #sleep(.5)
    if (dist(tool,end)<1):
      go = 0;
    
    
  return a,ang
    
def iteration(a, ang, tip_points, store_points):
  global ax
  
  # Clear the screen
  ax.clear()

  # Draw workspace
  plot_cuboid((19.24,0,15.24),(30.48,30.48,30.48))

  # Draw paper
  ax.plot(x,y,z, color='k')

  # Draw robot arm
  a.plot3D(ang)

  #print "Angles: ",ang
  #print("Tool tip position:")
  #print(a.getTool(ang)[0:3])
  if (store_points):
    tip_points[0].append(a.getTool(ang)[0])
    tip_points[1].append(a.getTool(ang)[1])
    tip_points[2].append(a.getTool(ang)[2])

  # Draw previous tool positions
  ax.plot_wireframe(tip_points[0], tip_points[1], tip_points[2], color='r')
  #draw()
  # Draw all buffered plots
  plt.draw()
  plt.pause(.001)
  #plt.iooff()
  show()
    
def main():
  global fig, ax,xend,yend,zend,x,y,z
  
  fig = gcf()
  ax = fig.gca(projection='3d')
  
  paper = PAPER  #4x4 rigid body transformation for paper position
  rotation = (paper[0:3,0:3])
  translation = paper[0:3,3:4]
  print(translation)
  point1 = rotation * np.matrix([[0],[0],[0]]) + translation
  point2 = rotation * np.matrix([[20.32],[0],[0]]) + translation
  point3 = rotation * np.matrix([[20.32],[27.94],[0]]) + translation
  point4 = rotation * np.matrix([[0],[27.94],[0]]) + translation
  print(point1)
  print(point2)
  print(point3)
  print(point4)

  basis1 = point2-point1
  basis1 = basis1/np.sqrt(np.multiply(basis1,basis1)).sum()
  #print('basis1')
  #print(basis1)
  basis2 = point4-point1
  basis2 = basis2/np.sqrt(np.multiply(basis2,basis2)).sum()

  convertpage = ConvertPage(point1, basis1, basis2)

  convertpage(0,0)

  x = asarray([float(point1[0]),float(point2[0]),float(point3[0]),float(point4[0]),float(point1[0])])
  y = asarray([float(point1[1]),float(point2[1]),float(point3[1]),float(point4[1]),float(point1[1])])
  z = asarray([float(point1[2]),float(point2[2]),float(point3[2]),float(point4[2]),float(point1[2])])
  print("paper")
  print(x)
  print(y)
  print(z)
  plt.show()
  
  a = Arm()
  
  xi = 30
  yi = 0
  zi = 35
  
  ang = a.angFromEnd(xi,yi,zi)

  
  tip_points = []
  tip_points.append([])
  tip_points.append([])
  tip_points.append([])

  strokes = PEN_STROKES
  store_points = False
  
  xend = 35
  yend = 0
  zend = 35

  while 1:
    

    # Get user input
    d = input("direction as list / angles as tuple?>")
    if type(d) == list:
      #a,ang = goToPoint(a,ang,d, tip_points, True)
      #end = asarray([[xend],[yend],[zend]])
      end = asarray([[d[0]],[d[1]],[d[2]]])
      a,ang = goToPoint(a,ang,end,tip_points,store_points)
      #Jt = a.getToolJac(ang)
      #ang = ang + dot(pinv(Jt)[:,:len(d)],d)
    elif type(d) == tuple:
      ang = d
      #iteration(a, ang, tip_points, store_points)
        
      #goToPoint(a,ang,asarray([[xend],[yend],[zend]]),tip_points,store_points)
    else:
      if (d == "reset"):
        store_points = False
        ang = [0,pi/4,pi/4]
        iteration(a, ang, tip_points, store_points)
      if (d == "clear" or "reset"):
        del tip_points[:][:]
        iteration(a, ang, tip_points, store_points)
      if (d == "draw"):
        for pts in strokes:
          print(pts)
          a,ang = goToPoint(a,ang,convertpage(pts[0][0]/10.0, pts[0][1]/10.0),tip_points, store_points)
          store_points = True
          sleep(1)
          
          a,ang = goToPoint(a,ang,convertpage(pts[1][0]/10.0, pts[1][1]/10.0),tip_points, store_points)
          
          sleep(1)
        store_points = False
  
  

  #a.plot3D(ang)
  #plt.draw()
  print(ang)
  show()
  
  
  
main()