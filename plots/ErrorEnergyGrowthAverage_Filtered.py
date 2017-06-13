import matplotlib.pyplot as plt
import numpy as np

k0 = 10
k1 = 10
tstart = [4,6,8,10,12]

ener = np.zeros((41,3,5))

for j in range(5):
    for m in range(3):
        if m == 0:
            infile = open("perturb_1024_h10_%0.2um/spec_kt" %tstart[j], "r")
        if m == 1:
            infile = open("perturb_1024_h10_%0.2ums/spec_kt" %tstart[j], "r")
        if m == 2:
            infile = open("perturb_1024_h10_%0.2us/spec_kt" %tstart[j], "r")

        enerlist = []
        for i in range(41):
            try:
                line = infile.readline()
                enerf = 0
                n = 0
                while n <= 512:
                    line = infile.readline()
                    if (n >= k0 and n <= k1):
                        words = line.split()
                        enerf += float(words[1])
                    n += 1
                enerlist.append(enerf)
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
plt.ylabel("Energy of Error at Scales $K \in [%d,%d]$" %(k0,k1))
plt.savefig("Error_Energy_Growth_Average_Filtered.pdf")
plt.close()
