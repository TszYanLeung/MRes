import matplotlib.pyplot as plt
import numpy as np

tstart = 9

ener = np.zeros((203,3))

for m in range(3):
    if m == 0:
        infile = open("perturb_1024_h10_%0.2ur_short/spec_t" %tstart, "r")
    if m == 1:
        infile = open("perturb_1024_h10_%0.2ul_short/spec_t" %tstart, "r")
    if m == 2:
        infile = open("perturb_1024_h10_%0.2us_short/spec_t" %tstart, "r")

    enerlist = []
    for i in range(203):
        try:
            line = infile.readline()
            words = line.split()
            enerlist.append(float(words[1]))
        except:
            pass
    ener[:len(enerlist),m] = np.asarray(enerlist)

    if len(enerlist) < 203:
        ener[len(enerlist):,m] = np.zeros(203-len(enerlist))

    if m == 0:
        plt.semilogy(np.arange(203)*0.0049444, ener[:,m], '-g')
    if m == 1:
        plt.semilogy(np.arange(203)*0.0049444, ener[:,m], '-b')
    if m == 2:
        plt.semilogy(np.arange(203)*0.0049444, ener[:,m], '-r')

    infile.close()

plt.xlabel("Lead Time")
plt.ylabel("Total Energy of Error")
plt.savefig("Error_Energy_Growth_%0.2u_Short.pdf" %tstart)
plt.close()
