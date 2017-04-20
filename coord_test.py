#Testing definition of mapping four corners of paper into a rigid body transform
#Worked for simple perfect case, need to implement with code given by Professor Revzen 
#  transform generation by least squares


#coordinate transform test
import numpy as np

from util import rigidFix

desired_result = np.matrix([[ 3./5, 4./5, 0.],[-4./5, 3./5, 0.],[0.,0.,1.]])


p1 = np.matrix([[0],[0],[0]]) 
p2 = np.matrix([[20.32],[0],[0]])
p3 = np.matrix([[20.32],[27.94],[0]]) 
p4 = np.matrix([[0],[27.94],[0]])

#r1 = np.matrix([[0],[0],[0]]) 
#r2 = np.matrix([[12.192],[-16.256],[0]])
#r3 = np.matrix([[34.544],[0.508],[0]]) 
#r4 = np.matrix([[22.352],[16.764],[0]])


r1 = np.matrix([[-.85],[50.226],[9.95]]) 
r2 = np.matrix([[-.22],[48.663],[22.217]])
r3 = np.matrix([[-15.73],[41.147],[23.97]]) 
r4 = np.matrix([[-17.57],[43.152],[10.96]])

'''
B = np.matmul(p1,r1.T)+np.matmul(p2,r2.T)+np.matmul(p3,r3.T)+np.matmul(p4,r4.T)
#B = np.matmul(p1,r1.T)+np.matmul(p2,r2.T)+np.matmul(p3,r3.T)

print(B)
U, s, V = np.linalg.svd(B, full_matrices=True)

M = np.diag(np.matrix([1,1,np.linalg.det(U)*np.linalg.det(V)]))

M = np.matrix([[1,0,0],[0,1,0],[0,0,np.linalg.det(U)*np.linalg.det(V)]])

print(M)
R = np.matmul(U,np.matmul(M,V))
print(desired_result)
R_i = np.linalg.inv(R)

print(R_i)
'''

P = p1

P = np.concatenate((P, p2), axis=1)
P = np.concatenate((P, p3), axis=1)
P = np.concatenate((P, p4), axis=1)



R = r1

R = np.concatenate((R, r2), axis=1)
R = np.concatenate((R, r3), axis=1)
R = np.concatenate((R, r4), axis=1)



#P = p1,p2,p3,p4
#R = r1,r2,r3,r4
#print(P)
#print(R)

#print(R.shape)

#print(desired_result)

print 

rot,ctr = rigidFix(P,R)