import numpy as np
import pandas as pd
from scipy.linalg import expm
import seaborn as sns
import matplotlib.pyplot as plt

T=2.5
delT = 1e-2
h = int(T/delT)

Q=np.array([[1e2, 0, 0],
              [0, 1, 0],
              [0, 0, 1]])

R=np.array([[19.9]])

x=np.zeros((3,h),dtype=np.float64)

u=np.zeros((1,h),dtype=np.float64)

P=np.zeros((3,3,h),dtype=np.float64)


a2 = 6.72e-3
m = 22.6e-3
g=9.81
L0=24.9*1e-3
L=0.520

x1eq = 4.5e-3
x2eq = 0
x3eq = np.sqrt((1/L0) * m * g * (2 * a2 * (1 + (x1eq / a2))**2))


k1 = L0 * x3eq / (a2 * (1 + x1eq / a2)**2)
k2 = L0 * x3eq**2 / (a2**2 * (1 + x1eq / a2)**3)


A=np.array([[0, 1, 0],
               [k2/m, 0, -k1/m],
               [0, 0, -R[0,0]/L]])

B = np.array([[0],
              [0],
              [1/L]])

Bt=np.transpose(B, axes=None)

At=np.transpose(A, axes=None)

R=np.array([[0.01]])

Ri= np.linalg.inv(R)

#resolver edo

H = np.block([[A, -B@Ri@Bt], [-Q, -At]])

for k in range(h-1,-1,-1):
    He = expm(-H*delT)
    XLambda = He@np.block([[np.identity(3)],[P[:,:,k]]])
    X = XLambda[0:3,0:3]
    Lambda = XLambda[3:6,0:3]
    Xi=np.linalg.inv(X)
    P[:,:,k-1]=Lambda@Xi

x[:,0]=[1e-3,0,0]

for k in range(h-1):
    u[:,k]=-Ri@Bt@P[:,:,k]@x[:,k]
    dxdt = A@x[:,k] + B@u[:,k]
    x[:,k+1] = x[:,k] + dxdt*delT


base = pd.DataFrame(data=np.transpose(x), columns=['serie_1', 'serie_2', 'serie_3'])
base.to_excel(r"C:\Users\DESKTOP\OneDrive\Documentos\Python\lq continuo.xlsx")

