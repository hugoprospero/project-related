import numpy as np
from scipy.linalg import expm

a=10

Q=np.array([[1, 0],
              [0, 1]])

R=np.array([[1]])

x=np.zeros((a,2),dtype=float)

u=np.zeros((a,1),dtype=float)

P=np.zeros((a,2,2),dtype=float)

A=np.array([[0.5, 1],
              [-0.1, 0.2]])

B = np.array([[0],
              [1]])

Bt=np.transpose(B, axes=None)

At=np.transpose(A, axes=None)

H= np.zeros((a,1,1),dtype=float)

Hi= np.zeros((a,1,1),dtype=float)



for k in range(a-1,0,-1):
    H[k] = R + Bt@P[k]@B
    Hi[k]=np.linalg.inv(H[k])
    P[k-1] = Q + At@(P[k]-P[k]@B@Hi[k]@Bt@P[k])@A

x[0]=[6,10]

for k in range(a-1):
    u[k]=-Hi[k]@Bt@P[k]@A@x[k]
    x[k+1]=A@x[k]+B@u[k]

print(x)
print(u)