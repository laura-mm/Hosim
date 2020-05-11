from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
import math
#plt.rc('text', usetex=True)
#plt.rcParams["text.latex.preamble"].append(r'\usepackage{xfrac}')
#plt.rc('xtick',labelsize=20)
#plt.rc('ytick',labelsize=20)
#plt.rc('axes', labelsize=20)
#plt.rc('axes', titlesize=20)
#plt.rc('legend', fontsize=20)
col = ['m', 'r', 'y', 'g', 'c', 'b']
colo = ['m.', 'r.', 'y.', 'g.', 'c.', 'b.']
lab = ['$\gamma = -\frac{1}{p-1}$', '$\gamma = 0$', '$\gamma = 1$']
mark = ['s', 'D', '^', 'v', '<', '>']
#lett = ['$(a)$', '$(b)$', '$(c)$', '$(d)$', '$(e)$', '$(f)$']

d = np.loadtxt("div23.txt", delimiter = ",")
d = d.reshape(len(d)/3, 3)
d = np.transpose(d)

#fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')

# Plot the surface
plt.plot_surface(d[0], d[1], d[2], color='b')

plt.figure(2)
plt.plot(d[1], d[2])

plt.show()

