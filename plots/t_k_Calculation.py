import matplotlib.pyplot as plt
import numpy as np
import math

tstart = [4,6,8,10,12]
run = "s"
resolution = 1.4
noctaves = int(math.log(512,resolution))

tdiff = np.zeros((len(tstart),noctaves))

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

    X = np.zeros(noctaves+2)
    Y = np.zeros(noctaves+1)
    for n in range(1,513):
        X[int(math.log(n,resolution))] += ener1[n]

    for i in range(noctaves+1):
        Y[i] = 0.5*(X[i]+X[i+1])
#        Y[i] = 0.5*(X[i])

    Z = np.zeros((41,noctaves+1))

    infile = open("perturb_1024_h10_%0.2u%s/spec_kt" %(tstart[m],run), "r")

    for i in range(41):
        line = infile.readline()
        ener2 = []
        n = 0
        while n <= 512:
            line = infile.readline()
            words = line.split()
            ener2.append(float(words[1]))
            n += 1
        ener2 = np.asarray(ener2)

        for n in range(1,513):
            Z[i,int(math.log(n,resolution))] += ener2[n]

    infile.close()

    t = np.zeros(noctaves+1)

    for j in range(noctaves+1):
        for i in range(41):
            if Z[i,j] < Y[j]:
                continue
            elif i > 0:
                t[j] = 0.2 * (i - (Z[i,j]-Y[j])/(Z[i,j]-Z[i-1,j]))
                break
            else:
                t[j] = 0
                break

    tdiff[m,:] = np.maximum(t[:-1] - t[1:],np.zeros(noctaves))
    #tdiffratio = np.divide(tdiff[m,:-1],tdiff[m,1:])

counter = np.zeros(noctaves)
for m in range(len(tstart)):
    for j in range(noctaves):
        if tdiff[m,j] != 0:
            counter[j] += 1

tdiffavg = np.sum(tdiff, axis=0)
for j in range(noctaves):
    if counter[j] != 0:
        tdiffavg[j] /= counter[j]
    else:
        tdiffavg[j] = 0

plt.semilogy(range(noctaves-1),tdiffavg[:-1],'-r')

def f(k):
    return 800*np.power(resolution**(-2.1),k)

plt.semilogy(range(8,14),f(range(8,14)),'-b')
plt.xlabel("$k$")
plt.ylabel("$t_k - t_{k+1}$ ($t_k = $ saturation time for $K \in [%.1f^k,%.1f^{k+1})$)" %(resolution,resolution))
plt.savefig("tdiff.pdf")
