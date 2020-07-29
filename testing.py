import numpy as np
import matplotlib.pyplot as plt

z = np.loadtxt('testing.txt', delimiter = ",")
z = z.reshape(len(z)//5, 5)
z = np.transpose(z)
"""
plt.figure(1)
plt.plot(z[0], z[1])
plt.plot(z[0], z[2])
plt.plot(z[0], z[3])
plt.ylim([-10, 10])
#plt.xlim([-1, 0])
plt.xlabel("z1")
plt.ylabel("sig2, sig3, M")
"""
plt.figure(2)
plt.plot(z[1], z[2])
plt.plot(z[1], z[3])
plt.plot(z[1], z[4])
#plt.xlim([0, 2])
plt.ylim([0, 5])
plt.xlabel("sig2")
plt.ylabel("sig3, M")

"""
#plt.plot(z[0], z[2])
#plt.axhline(y = 0, linestyle = '--', color = 'r')
#plt.ylim([-10, 10])
#plt.xlim([-10, 10])
"""

plt.show()



#file << zstart << "," << sig2 << "," << sig3(z1, help, sig2) << "," << sigtot(z1, help) << "," << M(zstart, help);
