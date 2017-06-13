import matplotlib.pyplot as plt
import numpy as np

infile = open("../../control_1024_h10_0x/spec_kt", "r")

for i in range(57*514):
    line = infile.readline()

k1 = np.arange(512)

line = infile.readline()
enst1a = []
n = 0
while n <= 511:
    line = infile.readline()
    words = line.split()
    enst1a.append(float(words[2]))
    n += 1
line = infile.readline()
enst1a = np.asarray(enst1a)

line = infile.readline()
enst1b = []
n = 0
while n <= 511:
    line = infile.readline()
    words = line.split()
    enst1b.append(float(words[2]))
    n += 1
line = infile.readline()
enst1b = np.asarray(enst1b)

line = infile.readline()
enst1c = []
n = 0
while n <= 511:
    line = infile.readline()
    words = line.split()
    enst1c.append(float(words[2]))
    n += 1
line = infile.readline()
enst1c = np.asarray(enst1c)

line = infile.readline()
enst1d = []
n = 0
while n <= 511:
    line = infile.readline()
    words = line.split()
    enst1d.append(float(words[2]))
    n += 1
enst1d = np.asarray(enst1d)

infile.close()

enst1 = 2*(enst1a + enst1b + enst1c + enst1d)*0.25

plt.loglog(k1, enst1, '-r')
plt.xlabel("Wavenumber")
plt.ylabel("Enstrophy Density")

infile = open("../spec_kt", "r")

for i in range(41):
    line = infile.readline()
    k2 = []
    enst2 = []
    n = 0
    while n <= 511:
        line = infile.readline()
        words = line.split()
        try:
            enst2.append(float(words[2]))
            k2.append(float(words[0]))
        except:
            pass
        n += 1
    line = infile.readline()

    if i % 5 == 0:
        plt.loglog(k2, enst2, '-b')

plt.savefig("Error_Enstrophy.pdf")
plt.close()

infile.close()
