#coordinate transform test
import numpy as np

desired_result = np.matrix([[ 3./5, 4./5, 0.],[-4./5, 3./5, 0.],[0.,0.,1.]])


p1 = np.matrix([[0],[0],[0]]) 
p2 = np.matrix([[20.32],[0],[0]])
p3 = np.matrix([[20.32],[27.94],[0]]) 
#p4 = np.matrix([[0],[27.94],[0]])

r1 = np.matrix([[0],[0],[0]]) 
r2 = np.matrix([[12.192],[-16.256],[0]])
r3 = np.matrix([[34.544],[0.508],[0]]) 
#r4 = np.matrix([[22.352],[16.764],[0]])

#B = np.matmul(p1,r1.T)+np.matmul(p2,r2.T)+np.matmul(p3,r3.T)+np.matmul(p4,r4.T)
B = np.matmul(p1,r1.T)+np.matmul(p2,r2.T)+np.matmul(p3,r3.T)

print(B)
U, s, V = np.linalg.svd(B, full_matrices=True)

M = np.diag(np.matrix([1,1,np.linalg.det(U)*np.linalg.det(V)]))

M = np.matrix([[1,0,0],[0,1,0],[0,0,np.linalg.det(U)*np.linalg.det(V)]])

print(M)
R = np.matmul(U,np.matmul(M,V))
print(desired_result)
R_i = np.linalg.inv(R)

print(R_i)