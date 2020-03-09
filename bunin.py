import numpy as np
import matplotlib.pyplot as plt

z = np.loadtxt('bunin2_-10.txt', delimiter = ",")
z = z.reshape(len(z)/3, 3)
z = np.transpose(z)

plt.figure(1)
plt.plot(z[0], z[1])
plt.plot(z[0], z[2])
plt.axhline(y = 0, linestyle = '--', color = 'r')

z = np.loadtxt('bunin2_-5.txt', delimiter = ",")
z = z.reshape(len(z)/3, 3)
z = np.transpose(z)

plt.figure(1)
plt.plot(z[0], z[1])
plt.plot(z[0], z[2])
plt.axhline(y = 0, linestyle = '--', color = 'r')

z = np.loadtxt('bunin2_0.txt', delimiter = ",")
z = z.reshape(len(z)/3, 3)
z = np.transpose(z)

plt.figure(1)
plt.plot(z[0], z[1])
plt.plot(z[0], z[2])
plt.axhline(y = 0, linestyle = '--', color = 'r')

z = np.loadtxt('bunin2_5.txt', delimiter = ",")
z = z.reshape(len(z)/3, 3)
z = np.transpose(z)

plt.figure(1)
plt.plot(z[0], z[1])
plt.plot(z[0], z[2])
plt.axhline(y = 0, linestyle = '--', color = 'r')

z = np.loadtxt('bunin2_10.txt', delimiter = ",")
z = z.reshape(len(z)/3, 3)
z = np.transpose(z)

plt.figure(1)
plt.plot(z[0], z[1])
plt.plot(z[0], z[2])
plt.axhline(y = 0, linestyle = '--', color = 'r')

plt.xlim([-5, 3])
plt.ylim([0, 5])

plt.show()
