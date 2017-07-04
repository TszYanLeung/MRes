import matplotlib.pyplot as plt
import numpy as np

tstart = 21

ener = np.zeros((512,31,4))

for m in range(4):
    if m == 0:
        infile = open("perturb_1024_h10_%0.2ur/spec_kt" %tstart, "r")
    if m == 1:
        infile = open("perturb_1024_h10_%0.2ul/spec_kt" %tstart, "r")
    if m == 2:
        infile = open("perturb_1024_h10_%0.2us/spec_kt" %tstart, "r")
    if m == 3:
        infile = open("perturb_1024_h10_%0.2us2/spec_kt" %tstart, "r")

    for i in range(31):
        enerlist = []
        line = infile.readline()
        line = infile.readline()
        for k in range(512):
            try:
                line = infile.readline()
                words = line.split()
                enerlist.append(words[1])
            except:
                pass
        ener[:len(enerlist),i,m] = np.asarray(enerlist)
        if len(enerlist) < 512:
            ener[len(enerlist):,i,m] = np.zeros(512-len(enerlist))

    growth = np.divide((ener[:,2:,m] - ener[:,1:-1,m])/0.3, ener[:,1:-1,m])
    inf = np.amin(growth[4:29], axis=0)
    sup = np.amax(growth[4:29], axis=0)
    infsum = np.zeros(30)
    supsum = np.zeros(30)
    for i in range(30):
        infsum[i] = np.sum(inf[:i])
        supsum[i] = np.sum(sup[:i])

    if m == 0:
        plt.plot(np.arange(29)*0.3+0.45, inf, '--g')
        plt.plot(np.arange(29)*0.3+0.45, sup, '--g')
        plt.plot(np.arange(30)*0.3+0.3, infsum, '-g')
        plt.plot(np.arange(30)*0.3+0.3, supsum, '-g')
    if m == 1:
        plt.plot(np.arange(29)*0.3+0.45, inf, '--b')
        plt.plot(np.arange(29)*0.3+0.45, sup, '--b')
        plt.plot(np.arange(30)*0.3+0.3, infsum, '-b')
        plt.plot(np.arange(30)*0.3+0.3, supsum, '-b')
    if m == 2:
        plt.plot(np.arange(29)*0.3+0.45, inf, '--r')
        plt.plot(np.arange(29)*0.3+0.45, sup, '--r')
        plt.plot(np.arange(30)*0.3+0.3, infsum, '-r')
        plt.plot(np.arange(30)*0.3+0.3, supsum, '-r')
    if m == 3:
        plt.plot(np.arange(29)*0.3+0.45, inf, '--m')
        plt.plot(np.arange(29)*0.3+0.45, sup, '--m')
        plt.plot(np.arange(30)*0.3+0.3, infsum, '-m')
        plt.plot(np.arange(30)*0.3+0.3, supsum, '-m')

    infile.close()

plt.plot(np.arange(30)*0.3+0.3, np.zeros(30), '-k')
plt.xlabel("Lead Time")
plt.ylabel("Growth Rate of Error (dashed) / Integral of Growth Rate (solid)")
plt.axis([0,9,-20,50])
plt.savefig("Growth_Rate_%0.2u.pdf" %tstart)
plt.close()
