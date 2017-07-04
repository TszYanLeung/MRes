import matplotlib.pyplot as plt
import numpy as np

tstart = 9

ener = np.zeros((512,203,3))

for m in range(3):
    if m == 0:
        infile = open("perturb_1024_h10_%0.2ur_short/spec_kt" %tstart, "r")
    if m == 1:
        infile = open("perturb_1024_h10_%0.2ul_short/spec_kt" %tstart, "r")
    if m == 2:
        infile = open("perturb_1024_h10_%0.2us_short/spec_kt" %tstart, "r")

    for i in range(203):
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

    growth = np.divide((ener[:,2:,m] - ener[:,1:-1,m])/0.0049444, ener[:,1:-1,m])
    inf = np.amin(growth, axis=0)
    sup = np.amax(growth, axis=0)
    infsum = np.zeros(202)
    supsum = np.zeros(202)
    for i in range(202):
        infsum[i] = np.sum(inf[:i])
        supsum[i] = np.sum(sup[:i])

    if m == 0:
        plt.plot(np.arange(201)*0.0049444+0.0074166, inf, '--g')
        plt.plot(np.arange(201)*0.0049444+0.0074166, sup, '--g')
        plt.plot(np.arange(202)*0.0049444+0.0049444, infsum, '-g')
        plt.plot(np.arange(202)*0.0049444+0.0049444, supsum, '-g')
    if m == 1:
        plt.plot(np.arange(201)*0.0049444+0.0074166, inf, '--b')
        plt.plot(np.arange(201)*0.0049444+0.0074166, sup, '--b')
        plt.plot(np.arange(202)*0.0049444+0.0049444, infsum, '-b')
        plt.plot(np.arange(202)*0.0049444+0.0049444, supsum, '-b')
    if m == 2:
        plt.plot(np.arange(201)*0.0049444+0.0074166, inf, '--r')
        plt.plot(np.arange(201)*0.0049444+0.0074166, sup, '--r')
        plt.plot(np.arange(202)*0.0049444+0.0049444, infsum, '-r')
        plt.plot(np.arange(202)*0.0049444+0.0049444, supsum, '-r')

    infile.close()

plt.plot(np.arange(202)*0.0049444+0.0049444, np.zeros(202), '-k')
plt.xlabel("Lead Time")
plt.ylabel("Growth Rate of Error (dashed) / Integral of Growth Rate (solid)")
#plt.axis([0,1,-20,100])
plt.savefig("Growth_Rate_%0.2u_Short.pdf" %tstart)
plt.close()
