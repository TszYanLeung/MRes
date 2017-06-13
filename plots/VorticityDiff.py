import matplotlib.pyplot as plt
import numpy as np

# Initialise the coordinates
x = np.arange(0, 1024)
y = np.arange(0, 1024)

for i in range(20,61,5):
    # Load the data and put it into the correct shape
    data1 = np.loadtxt("../qq%0.4u.dat" %i, skiprows=1)
    data2 = np.loadtxt("../../control_1024_h10_0x/qq%0.4u.dat" %i, skiprows=1)
    z1 = np.transpose(np.reshape(data1, (1024, 1024)))
    z2 = np.transpose(np.reshape(data2, (1024, 1024)))
    z = z1 - z2

    # Plotting the vorticity contours
    plt.contourf(x, y, z)
    plt.xlabel("x")
    plt.ylabel("y")
    plt.savefig("Error_Vorticity_%.2f.pdf" %(0.2*i))
    plt.close()
