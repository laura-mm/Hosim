import numpy as np
import matplotlib.pyplot as plt

z = np.loadtxt('testing.txt', delimiter = ",")
z = z.reshape(len(z)/3, 3)
z = np.transpose(z)

plt.figure(1)
plt.plot(z[0], z[1])
#plt.plot(z[0], z[2])
plt.axhline(y = 0, linestyle = '--', color = 'r')
plt.ylim([-10, 10])
plt.show()
