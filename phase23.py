from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
import math
data = np.loadtxt("crit23.txt", delimiter = ",")
data = data.reshape(1001, 1001, 4) # mu2 mu3 (mu2, mu3, sig2, sig3)
data = np.transpose(data, (2, 0, 1)) # (mu2, mu3, sig2, sig3), mu2, mu3

mu2 = data[0]
mu3 = data[1]
sig2 = data[2]
sig3 = data[3]

fig = plt.figure()
ax = fig.gca(projection = '3d')

surf1 = ax.plot_surface(mu2, mu3, sig2)
surf2 = ax.plot_surface(mu2, mu3, sig3)
plt.show()

"""
'''
======================
3D surface (color map)
======================

Demonstrates plotting a 3D surface colored with the coolwarm color map.
The surface is made opaque by using antialiased=False.

Also demonstrates using the LinearLocator and custom formatting for the
z axis tick labels.
'''

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np


fig = plt.figure()
ax = fig.gca(projection='3d')

# Make data.
X = np.arange(-5, 5, 0.25)
Y = np.arange(-5, 5, 0.25)
X, Y = np.meshgrid(X, Y)
R = np.sqrt(X**2 + Y**2)
Z = np.sin(R)

# Plot the surface.
surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)

# Customize the z axis.
ax.set_zlim(-1.01, 1.01)
ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

# Add a color bar which maps values to colors.
fig.colorbar(surf, shrink=0.5, aspect=5)

plt.show()
"""
