# heat maps and phase diagrams

import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import math
import matplotlib.colors as colors
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
#plt.rc('text', usetex=True)
#plt.rcParams["text.latex.preamble"].append(r'\usepackage{xfrac}')
lett = ['$(a)$', '$(b)$', '$(c)$']

plt.rc('xtick',labelsize=20)
plt.rc('ytick',labelsize=20)
plt.rc('axes', labelsize=20)
plt.rc('axes', titlesize=20)

grid = 100 # for heat map
gridl = 1000 # lines on heat map
gridg = 200 # 3d plots

# could also look at interpolation - to make it more smooth and less pixelly

"""
# fp solution plot acrit only
################################################

ast = 3 #start for a
gst = 3 # start for gamma

plt.figure(1)
ag = np.linspace((gst*(2.0/gridg)) - 1, 1, gridg + 1 - gst)
a = np.linspace(ast*(2.0/gridg), 2, gridg + 1 - ast)
ag, a = np.meshgrid(ag, a)
acrit = np.loadtxt('acrit_0.txt', delimiter=",")
acrit = acrit.reshape(gridg + 1 - gst, gridg + 1 - ast)
acrit = np.transpose(acrit)
plt.pcolormesh(ag, a, acrit, norm = colors.LogNorm(vmin = acrit.min(), vmax = acrit.max()), cmap = cm.gist_ncar)
plt.axhline(y = 1, linestyle = '--', color = 'r')
plt.colorbar()
plt.xlabel(r'$\gamma$')
plt.ylabel(r'$a$')
plt.xlim([-1, 1])
plt.ylim([0, 2])
plt.subplots_adjust(left = 0.1, right = 1.0, top = 0.96, bottom = 0.1)
plt.savefig('acrit.pdf')

###################################
"""
#heat maps
##############################################

g = np.linspace(-1, 1, grid + 1)
x = np.linspace(-1, 3, grid + 1)
s = np.power(10, x)

#la5 = np.loadtxt('phline_5_0.txt', delimiter=",")
#la2 = np.loadtxt('phline_20_0.txt', delimiter=",")

#gl = np.linspace((2.0/gridl) - 1, 1, gridl)
#l1 = math.sqrt(2)/(1+gl);


# order 2 heatmaps, 2 values for mu
####################################################

plt.subplots(1, 2, figsize=(12,6))

plt.subplot(1, 2, 1) # mu 0
A0 = np.loadtxt('hoA0_lt.txt', delimiter=",")
A0 = A0.reshape(grid + 1, grid + 1, 3)
A0 = np.transpose(A0, (1, 0, 2))
plt.imshow(A0, extent=[-1, 1, 0.1, 1000], origin = 'lower')
plt.yscale('log')
#plt.semilogy(gl, l1, 'k')
#plt.xlabel(r'$\gamma$')
#plt.ylabel(r'$\sigma$')
plt.text(0.8, 500, lett[0], fontsize=20)
#plt.title(r'$a = 2$')

"""
plt.subplot(1, 2, 2) # p5
p5 = np.loadtxt('p5.txt', delimiter=",")
p5 = p5.reshape(grid + 1, grid + 1, 3)
p5 = np.transpose(p5, (1, 0, 2))
plt.imshow(p5, extent=[-1, 1, 0.1, 1000], origin = 'lower')
plt.yscale('log')
plt.semilogy(gl, la5, 'k')
plt.xlabel(r'$\gamma$')
plt.ylim([0.1, 1000])
plt.text(0.8, 500, lett[1], fontsize=20)
plt.title(r'$a = 0.5$')
"""

plt.subplots_adjust(left = 0.07, right = 0.98, wspace = 0.17)
plt.savefig('hoA.pdf')



###############################
# end of heat maps


plt.show()




