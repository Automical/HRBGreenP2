import numpy as np
from scipy.linalg import expm as expM
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

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


class Arm( object ):
  """
  class Arm
  
  Represents a series manipulator made of several segments.
  Each segment is graphically represented by a wireframe model
  """
  def __init__(self):
    # link lengths
    self.ll = asarray([10,10,10])
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
    #print('what')
    #print('tw')
    #print(tw)
    #print('ang')
    #print(ang)
    #print('selftw')
    #print(self.tw)
    A = [identity(4)]
    #print('calling at')
    for twi in tw:
      M = seToSE(twi)
      #print('twi')
      #print(twi)
      #print('M')
      #print(M)
      A.append(dot(A[-1],M))
    return A
    
  def getTool( self, ang ):
    """
    Get "tool tip" position in world coordinates
    """
    # Get the rigid transformation for the last segment of the arm
    M = self.at(ang)[-1]
    #print('M')
    #print(M)
    return dot(M, self.tool)
  
  def getToolJac( self, ang ):
    """
    Get "tool tip" Jacobian by numerical approximation
    
    NOTE: implementation is a placeholder. This method is overwritten
    dynamically by __init__() to point to a jacobian_cdas() function
    """
    raise RuntimeError("uninitialized method called")
    
  def plotIJ( self, ang, axI=0, axJ=1 ):
    """
    Display the specified axes of the arm at the specified set of angles
    """
    A = self.at(ang)
    print("A")
    print(A)
    print("geom")
    print(self.geom)
    print("zip")
    print(zip(A,self.geom))
    for a,g in zip(A, self.geom):
      ng = dot(a,g)
      plot( ng[axI,:], ng[axJ,:], '.-' )
    tp = dot(a, self.tool)
    plot( tp[axI], tp[axJ], 'hk' )
    plot( tp[axI], tp[axJ], '.y' )
    

  def plot3D( self, ang ):
    global ax
    """
    Plot arm in 3 views
    """
    '''
    #Axes3D.plot_wireframe(X, Y, Z)
    A = self.at(ang)
    print("A")
    print(A)
    print("shape A")
    print(shape(A))
    print("ang")
    print(ang)
    #print("zip")
    #print(zip(A,self.geom))
    #print("starting loop")
    print('tw')
    print(self.tw)
    print('get tool')
    t=self.getTool( ang );
    print(t)
    '''
    '''
    for a,g in zip(A, self.geom):
      print('a')
      print(a)
      print('g')
      print(g)
      ng = dot(a,g)
      print('ng0')
      print(ng[0,:])
      print('ng1')
      print(ng[1,:])
      print('ng2')
      print(ng[2,:])
      #print(ng)
      #print(shape(ng))
      Axes3D.plot(self.ax,ng[0,:],ng[1,:],ng[2,:])
      #plot( ng[axI,:], ng[axJ,:], '.-' )
      '''
    x = []
    y = []
    z = []
    #print('A')
    A = self.at(ang)
    #print("zip")
    #print(zip(A,self.geom))
    for a,g in zip(A, self.geom):
      #print('loop start')
      ng = dot(a,g)
      '''
      print(ng)
      print(shape(ng))
      print(shape(ng[0]))
      print(shape(ng[1]))
      print(shape(ng[2]))
      '''
      
      print('ng test')
      print('ng')
      print(ng)
      print('ng0')
      print(ng[0])
      print('ng1')
      print(ng[1])
      print('ng2')
      print(ng[2])
      print('x')
      print(x)
      print('y')
      print(y)
      print('z')
      print(z)
      x.append(float(ng[0]))
      y.append(float(ng[1]))
      z.append(float(ng[2]))
    #self.ax.plot( ng[0,:], ng[1,:], ng[2,:], label='parametric curve' )
    print(x)
    print(y)
    print(z)
    tp = self.getTool(ang)
    x.append(float(tp[0]))
    y.append(float(tp[1]))
    z.append(float(tp[2]))
    print(x)
    print(y)
    print(z)
    ax.plot( x, y, z,  zdir='z' )
    ax.plot( x, y, z, 'hk', zdir='z' )
    #ax.plot(tp[0],tp[1],tp[2],'hr',zdir='z')
    ax.set_xlim(-30,30)
    ax.set_xlabel('x')
    ax.set_ylim(-30,30)
    ax.set_ylabel('y')
    ax.set_zlim(-30,30)
    ax.set_zlabel('z')
    #tp = dot(a, self.tool)
    #tp = self.getTool(ang)
    tx = asarray([0,(tp[0])])
    ty = asarray([0,(tp[1])])
    tz = asarray([0,(tp[2])])
    #print(tp)
    #print(tx)
    #print(ty)
    #print(tz)
    #ax.plot( [x[-1],tp[0]], [y[-1],tp[1], [z[-1],tp[2]], 'hk' )
    #[x[-1],ax.pl]ot[y[-1],( tp[0], tp[1], tp[2], '.y' )
    ax.plot(tx,ty,tz,'hr',zdir='z')
    #print(A)
    #tp = dot(a, self.tool)
    #print("tool")
    #print(tp)
    
    '''
    A = self.at(ang)
    print("A")
    print(A)
    
    ax = [-20,20,-20,20]
    subplot(2,2,1)
    self.plotIJ(ang,0,1)
    axis('equal')
    axis(ax)
    grid(1)
    xlabel('X'); ylabel('Y')
    subplot(2,2,2)
    self.plotIJ(ang,2,1)
    axis('equal')
    axis(ax)
    grid(1)
    xlabel('Z'); ylabel('Y')
    subplot(2,2,3)
    self.plotIJ(ang,0,2)
    axis('equal')
    axis(ax)
    grid(1)
    xlabel('X'); ylabel('Z')
    '''

class ConvertPage(object):
  def __init__(self, o, b1, b2):
    self.origin = o
    self.basis1 = b1
    self.basis2 = b2
        
  def __call__(self, x, y):
    return self.origin + self.basis1 * x + self.basis2 * y


def main():
  global fig, ax
  """
  Run an example of a robot arm
  
  This can be steered via inverse Jacobian, or positioned.
  """
  fig = gcf()
  ax = fig.gca(projection='3d')

  paper = np.matrix([[ 3./5, 4./5, 0., 50.],[-4./5, 3./5, 0., 0.],[0.,0.,1.,0.],[0.,0.,0.,1.]])  #4x4 rigid body transformation for paper position
  rotation = paper[0:3,0:3]
  translation = paper[0:3,3:4]/10.0
  print(translation)
  point1 = rotation * np.matrix([[0],[0],[0]]) + translation
  point2 = rotation * np.matrix([[20.32],[0],[0]]) + translation
  point3 = rotation * np.matrix([[20.32],[27.94],[0]]) + translation
  point4 = rotation * np.matrix([[0],[27.94],[0]]) + translation

  basis1 = point2-point1
  basis1 = basis1/np.sqrt(np.multiply(basis1,basis1)).sum()
  #print('basis1')
  #print(basis1)
  basis2 = point3-point1
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
  #f = gcf()
  ang = [0,pi/4,pi/4]
  while 1:
    #a.fig.set(visible=0)
    #clf()
    ax.clear()
    ax.plot_wireframe(x,y,z)
    a.plot3D(ang)
    
    plt.draw()
    #a.fig.set(visible=1)
    #draw()
    print "Angles: ",ang
    d = input("direction as list / angles as tuple?>")
    if type(d) == list:
      Jt = a.getToolJac(ang)
      ang = ang + dot(pinv(Jt)[:,:len(d)],d)
    else:
      ang = d
  
main()
