import numpy as np
from scipy.linalg import expm as expM
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from time import sleep
from itertools import product, combinations

from math import *

def dist(x1,y1,x2,y2):
  return sqrt((x2-x1)**2+(y2-y1)**2);

class Arm( object ):

  def __init__(self):
    # link lengths in cm
    self.l1 = 28
    self.l2 = 35
    self.l3 = 20
    self.l4 = 15
    self.l5 = 35
    self.l6 = 10
    self.l7 = 10
    self.l8 = 5
    
  def angFromEnd(self,x,y):
    ang = [90,90,90,90];
    
    
    d1 = dist(self.l6,0,x,y)
    print("d1 = %f" % d1)
    beta = acos((self.l2**2+d1**2-self.l5**2)/(2*self.l2*d1))
    print("beta = %f" % degrees(beta))
    alpha = atan2(y,x-self.l6)
    print("alpha = %f" % degrees(alpha))
    theta2 = alpha + beta;
    print("theta2 = %f" % degrees(theta2))

    
    ang[1] = theta2
    
    xc = self.l2*cos(theta2) + self.l6
    yc = self.l2*sin(theta2)
    
    print("xc = %f" % xc)
    print("yc = %f" % yc)
    
    theta4 = -atan2(yc-y,x-xc)
    ang[3]=theta4;
    
    #theta4 = acos((x-self.l2*cos(theta2))/(self.l5))
    print("theta4 = %f" % degrees(theta4))
    #xb = x-(self.l4+self.l5)*cos(theta4)
    #yb = y-(self.l4+self.l5)*sin(theta4)
    xb = self.l6 + self.l2*cos(theta2) + self.l4*cos(theta4)
    yb = self.l2*sin(theta2) + self.l4*sin(theta4)
    print("xb = %f" % xb)
    print("yb = %f" % yb)
    d2 = dist(-self.l7,0,xb,yb)
    print("d2 = %f" % d2)
    gamma = acos((self.l1**2+d2**2-self.l3**2)/(2*self.l1*d2))
    print("gamma = %f" % degrees(gamma))
    delta = atan2(yb,xb+self.l7)
    print("delta = %f" % degrees(delta))
    theta1 = gamma+delta
    
    if (theta1<(pi/2)):
      theta1 = (pi/2-theta1) + pi/2
    
    print("theta1= %f" % degrees(theta1))
    ang[0] = theta1
    
    theta3 = acos((self.l1**2+self.l3**2-d2**2)/(2*self.l1*self.l3))
    ang[2] = theta3
    return ang
    
  def plot3D(self,ang):
    global ax
    xc = self.l2*cos(ang[1]) + self.l6
    yc = self.l2*sin(ang[1])
    xe = self.l5*cos(ang[3])+xc
    ye = self.l5*sin(ang[3])+yc
    xa = -self.l7+self.l1*cos(ang[0])
    ya = self.l1*sin(ang[0])
    xb = self.l3*cos(ang[2])+xa
    yb = self.l3*sin(ang[2])+ya
    
    x = [-self.l7, xa]
    y = [0,ya]
    z = [0,0]
    ax.plot(x,y,z,label='l1')
    x = [xa, xb]
    y = [ya, yb]
    ax.plot(x,y,z,label='l3')
    x = [xb,xc]
    y = [yb,yc]
    ax.plot(x,y,z,label='l4')
    x = [self.l6, xc]
    y = [0,yc]
    ax.plot(x,y,z,label='l2')
    x = [xc,xe]
    y= [yc,ye]
    ax.plot(x,y,z,label='l5')
    
    
    
def main():
  global fig, ax,x,y,z
  
  fig = gcf()
  ax = fig.gca(projection='3d')
  
  
  a = Arm()
  
  x = 40
  y = 20
  
  ang = a.angFromEnd(x,y)
  a.plot3D(ang)
  plt.draw()
  print(ang)
  show()
  
  
  
main()