import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

alpha = 2.0

k1 = 413
k2 = -47

A = np.zeros((1025,1025))
for i in range(1025):
    for j in range(1025):
        l1 = i-512
        l2 = j-512
        absk = (k1**2+k2**2)**0.5
        absl = (l1**2+l2**2)**0.5
        abskl = ((k1-l1)**2+(k2-l2)**2)**0.5
        A[i,j] = -(absk**(-alpha))*((abskl**alpha)-(absl**alpha))*(k1*l2-k2*l1)

x = np.arange(-512,513)
y = np.arange(-512,513)
plt.contourf(x,y,np.transpose(A))
plt.colorbar()
plt.savefig("Interaction_Coeff.pdf")
plt.close()
