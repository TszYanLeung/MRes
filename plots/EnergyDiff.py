import matplotlib.pyplot as plt
import numpy as np

infile = open("../../control_1024_h10_0x/spec_kt", "r")

for i in range(57*514):
    line = infile.readline()

k1 = np.arange(512)

line = infile.readline()
ener1a = []
n = 0
while n <= 511:
    line = infile.readline()
    words = line.split()
    ener1a.append(float(words[1]))
    n += 1
line = infile.readline()
ener1a = np.asarray(ener1a)

line = infile.readline()
ener1b = []
n = 0
while n <= 511:
    line = infile.readline()
    words = line.split()
    ener1b.append(float(words[1]))
    n += 1
line = infile.readline()
ener1b = np.asarray(ener1b)

line = infile.readline()
ener1c = []
n = 0
while n <= 511:
    line = infile.readline()
    words = line.split()
    ener1c.append(float(words[1]))
    n += 1
line = infile.readline()
ener1c = np.asarray(ener1c)

line = infile.readline()
ener1d = []
n = 0
while n <= 511:
    line = infile.readline()
    words = line.split()
    ener1d.append(float(words[1]))
    n += 1
ener1d = np.asarray(ener1d)

infile.close()

ener1 = 2*(ener1a + ener1b + ener1c + ener1d)*0.25

plt.loglog(k1, ener1, '-r')
plt.xlabel("Wavenumber")
plt.ylabel("Energy Density")

infile = open("../spec_kt", "r")

for i in range(41):
    line = infile.readline()
    k2 = []
    ener2 = []
    n = 0
    while n <= 511:
        line = infile.readline()
        words = line.split()
        try:
            ener2.append(float(words[1]))
            k2.append(float(words[0]))
        except:
            pass
        n += 1
    line = infile.readline()

    if i % 5 == 0:
        plt.loglog(k2, ener2, '-b')

plt.savefig("Error_Energy.pdf")
plt.close()

infile.close()
