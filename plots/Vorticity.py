import matplotlib.pyplot as plt
import numpy as np

# Initialise the coordinates
x = np.arange(0, 1024)
y = np.arange(0, 1024)

for i in range(0,101,5):
    # Load the data and put it into the correct shape
    data = np.loadtxt("../qq%0.4u.dat" %i, skiprows=1)
    z = np.transpose(np.reshape(data, (1024, 1024)))

    # Plotting the vorticity contours
    plt.contourf(x, y, z)
    plt.xlabel("x")
    plt.ylabel("y")
    plt.savefig("Vorticity_%.2f.pdf" %(0.2*i))
    plt.close()
