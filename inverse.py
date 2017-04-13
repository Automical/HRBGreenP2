import numpy as np
from scipy.linalg import expm as expM
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from time import sleep
from itertools import product, combinations

from math import *

class Arm( object ):

  def __init__(self):
    # link lengths in cm
    self.l1 = 10
    self.l2 = 10
    self.l3 = 10
    self.l4 = 5
    self.l5 = 10
    self.l6 = 8
    self.l7 = 8
    self.l8 = 5
    
  def angFromEnd(self,x,y):
    ang = [90,90];
    
    d1 = sqrt((x-self.l6)**2+y**2)
    print("d1 = %f\n" % d1)
    beta = acos((self.l2**2+d1**2-self.l5**2)/(2*self.l5*d1))
    print("beta = %f\n" % beta)
    alpha = atan2(y,x+self.l6)
    print("alpha = %f\n" % alpha)
    theta2 = alpha + beta;
    print("theta2 = %f\n" % theta2)

    
    ang[1] = theta2
    
    theta4 = acos((x-self.l2*cos(theta2))/(self.l5))
    
    xb = x-(self.l4+self.l5)*cos(theta4)
    yb = y-(self.l4+self.l5)*sin(theta4)
    
    d2 = sqrt((xb+self.l7)**2+self.l7**2)
    
    gamma = acos((self.l1**2+d2**2-self.l3**2)/(2*self.l1*d2))
    
    delta = atan2(yb,xb+self.l7)
    
    theta1 = gamma+delta
    
    ang[0] = theta1
    
    return ang
    
    
def main():
  
  a = Arm()
  
  ang = a.angFromEnd(15,10)
  
  print(ang)
  
  
main()