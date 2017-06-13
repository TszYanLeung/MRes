import matplotlib.pyplot as plt
import numpy as np

tstart = [4,6,8,10,12]

ener = np.zeros((41,3,5))

for j in range(5):
    for m in range(3):
        if m == 0:
            infile = open("perturb_1024_h10_%0.2um/spec_t" %tstart[j], "r")
        if m == 1:
            infile = open("perturb_1024_h10_%0.2ums/spec_t" %tstart[j], "r")
        if m == 2:
            infile = open("perturb_1024_h10_%0.2us/spec_t" %tstart[j], "r")

        enerlist = []
        for i in range(41):
            try:
                line = infile.readline()
                words = line.split()
                enerlist.append(float(words[1]))
            except:
                pass
        ener[:len(enerlist),m,j] = np.asarray(enerlist)

        if len(enerlist) < 41:
            ener[len(enerlist):,m,j] = np.zeros(41-len(enerlist))

        infile.close()

eneravg = 0.2*np.sum(ener, axis=2)

for m in range(3):
    if m == 0:
        plt.semilogy(np.arange(41)*0.2, eneravg[:,m], '-g')
    if m == 1:
        plt.semilogy(np.arange(41)*0.2, eneravg[:,m], '-b')
    if m == 2:
        plt.semilogy(np.arange(41)*0.2, eneravg[:,m], '-r')

plt.xlabel("Lead Time")
plt.ylabel("Total Energy of Error")
plt.savefig("Error_Energy_Growth_Average.pdf")
plt.close()
