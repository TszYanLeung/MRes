import matplotlib.pyplot as plt
import numpy as np

k0 = 5
tstart = 12

ener = np.zeros((41,3))

for m in range(3):
    if m == 0:
        infile = open("perturb_1024_h10_%0.2um/spec_kt" %tstart, "r")
    if m == 1:
        infile = open("perturb_1024_h10_%0.2ums/spec_kt" %tstart, "r")
    if m == 2:
        infile = open("perturb_1024_h10_%0.2us/spec_kt" %tstart, "r")

    enerlist = []
    for i in range(41):
        try:
            line = infile.readline()
            enerf = 0
            n = 0
            while n <= 512:
                line = infile.readline()
                if n >= k0:
                    words = line.split()
                    enerf += float(words[1])
                n += 1
            enerlist.append(enerf)
        except:
            pass
    ener[:len(enerlist),m] = np.asarray(enerlist)

    if len(enerlist) < 41:
        ener[len(enerlist):,m] = np.zeros(41-len(enerlist))

    if m == 0:
        plt.semilogy(np.arange(41)*0.2, ener[:,m], '-g')
    if m == 1:
        plt.semilogy(np.arange(41)*0.2, ener[:,m], '-b')
    if m == 2:
        plt.semilogy(np.arange(41)*0.2, ener[:,m], '-r')

    infile.close()

plt.xlabel("Lead Time")
plt.ylabel("Energy of Error at Scales $K > %d$" %k0)
plt.savefig("Error_Energy_Growth_Filtered_%0.2u.pdf" %tstart)
plt.close()
