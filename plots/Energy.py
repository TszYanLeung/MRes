import matplotlib.pyplot as plt
import numpy as np

def f(k):
    return 0.5*np.power(k, -8./3.)

def g(k):
    return 0.05*np.power(k, -2.)

infile = open("../spec_kt", "r")

for i in range(101):
    line = infile.readline()
    k = []
    qsp = []
    n = 0
    while n <= 511:
        line = infile.readline()
        words = line.split()
        try:
            qsp.append(float(words[1]))
            k.append(float(words[0]))
        except:
            pass
        n += 1
    line = infile.readline()

    if i % 5 == 0:
        plt.loglog(k, qsp, '-r')
        plt.loglog(k[10:100], f(k[10:100]), '-k')
        plt.loglog(k[100:300], g(k[100:300]), '-k')
        plt.xlabel("Wavenumber")
        plt.ylabel("Energy Spectral Density")
        plt.savefig("Energy_%.2f.pdf" %(0.2*i))
        plt.close()

infile.close()
