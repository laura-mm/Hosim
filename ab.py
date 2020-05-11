# this used to plot functions from numerical solution

import numpy as np
import matplotlib.pyplot as plt
N = 200000
x = np.loadtxt('abs.txt', delimiter=",")
x = x.reshape(len(x)/3, 3)
x = np.transpose(x)

plt.figure(1)
plt.plot(x[0], x[1])
plt.plot(x[0], x[2])
plt.ylim([-10, 10])
plt.axhline(y = 0, linestyle = '--', color = 'r')
plt.axvline(x = -0.839924, linestyle = '--', color = 'k')
plt.show()

