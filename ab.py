# this used to plot functions from numerical solution

import numpy as np
import matplotlib.pyplot as plt
N = 200000
x = np.loadtxt('zz.txt', delimiter=",")
y = np.loadtxt('abs.txt', delimiter=",")
y = y.reshape(len(y)/2, 2)
y = np.transpose(y)

plt.figure(1)
plt.plot(x, y[0])
#plt.scatter(x, y[0])
plt.plot(x, y[1])
plt.ylim([-10, 10])
plt.axhline(y = 0, linestyle = '--', color = 'r')
plt.show()

