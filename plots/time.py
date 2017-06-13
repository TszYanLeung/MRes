import matplotlib.pyplot as plt
import numpy as np
import math

tstart = [4,6,8,10,12]
run = "s"

t = np.zeros((len(tstart),512))

for m in range(len(tstart)):
    infile = open("control_1024_h10_0x/spec_kt", "r")

    for i in range((37+5*tstart[m])*514):
        line = infile.readline()

    line = infile.readline()
    ener1a = []
    n = 0
    while n <= 512:
        line = infile.readline()
        words = line.split()
        ener1a.append(float(words[1]))
        n += 1
    ener1a = np.asarray(ener1a)

    line = infile.readline()
    ener1b = []
    n = 0
    while n <= 512:
        line = infile.readline()
        words = line.split()
        ener1b.append(float(words[1]))
        n += 1
    ener1b = np.asarray(ener1b)

    line = infile.readline()
    ener1c = []
    n = 0
    while n <= 512:
        line = infile.readline()
        words = line.split()
        ener1c.append(float(words[1]))
        n += 1
    ener1c = np.asarray(ener1c)

    line = infile.readline()
    ener1d = []
    n = 0
    while n <= 512:
        line = infile.readline()
        words = line.split()
        ener1d.append(float(words[1]))
        n += 1
    ener1d = np.asarray(ener1d)

    infile.close()

    ener1 = 2*(ener1a + ener1b + ener1c + ener1d)*0.25

    Y = 0.55*ener1[:-1]

    Z = np.zeros((41,512))

    infile = open("perturb_1024_h10_%0.2u%s/spec_kt" %(tstart[m],run), "r")

    for i in range(41):
        line = infile.readline()
        ener2 = []
        n = 0
        while n <= 511:
            line = infile.readline()
            words = line.split()
            ener2.append(float(words[1]))
            n += 1
        line = infile.readline()
        ener2 = np.asarray(ener2)
        Z[i,:] = ener2

    infile.close()

    for j in range(512):
        for i in range(41):
            if Z[i,j] < Y[j]:
                continue
            elif i > 0:
                t[m,j] = 0.2 * (i - (Z[i,j]-Y[j])/(Z[i,j]-Z[i-1,j]))
                break
            else:
                t[m,j] = 0
                break

counter = np.zeros(512)
for m in range(len(tstart)):
    for j in range(512):
        if t[m,j] != 0:
            counter[j] += 1

tavg = np.sum(t, axis=0)
for j in range(512):
    if counter[j] != 0:
        tavg[j] /= counter[j]
    else:
        tavg[j] = 0

plt.loglog(range(512),tavg,'-r')

def f(k):
    return 90*np.power(k, -1.25)

plt.semilogy(range(6,201),f(range(6,201)),'-b')
plt.xlabel("Wavenumber")
plt.ylabel("Saturation Time")
plt.savefig("time.pdf")
