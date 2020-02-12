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

p2m0 = np.loadtxt('phline_2_0.txt', delimiter=",")
p2m0 = p2m0.reshape(int(len(p2m0)/2), 2)
p2m0 = np.transpose(p2m0)

p2m2 = np.loadtxt('phline_2_-20.txt', delimiter=",")
p2m2 = p2m2.reshape(int(len(p2m2)/2), 2)
p2m2 = np.transpose(p2m2)

#la2 = np.loadtxt('phline_20_0.txt', delimiter=",")

#gl = np.linspace((2.0/gridl) - 1, 1, gridl)
#l1 = math.sqrt(2)/(1+gl);


# order 2 heatmaps, 2 values for mu
####################################################

plt.subplots(1, 2, figsize=(12,6))

plt.subplot(1, 2, 1) # mu 0
A0 = np.loadtxt('ho20_lt.txt', delimiter=",")
A0 = A0.reshape(grid + 1, grid + 1, 3)
A0 = np.transpose(A0, (1, 0, 2))
#plt.imshow(A0, extent=[-1, 1, 0.1, 1000], origin = 'lower')
plt.imshow(A0, extent=[-1, 1, 1, 2], origin = 'lower')
plt.yscale('log')
#plt.semilogy(p2m0[0], p2m0[1], 'k')
#plt.xlabel(r'$\gamma$')
#plt.ylabel(r'$\sigma$')
#plt.ylim([0.1, 1000])
plt.text(0.8, 500, lett[0], fontsize=20)
#plt.title(r'$a = 2$')



plt.subplot(1, 2, 2) # mu -2
A2 = np.loadtxt('ho2-_lt.txt', delimiter=",")
A2 = A2.reshape(grid + 1, grid + 1, 3)
A2 = np.transpose(A2, (1, 0, 2))
plt.imshow(A2, extent=[-1, 1, 0.1, 1000], origin = 'lower')
#plt.imshow(A2, extent=[-1, 1, 1, 1000], origin = 'lower')
plt.yscale('log')
plt.semilogy(p2m2[0], p2m2[1], 'k')
#plt.plot(p2m2[0], np.log10(p2m2[1]), 'k')
#plt.xlabel(r'$\gamma$')
plt.ylim([0.1, 1000])
plt.text(0.8, 500, lett[1], fontsize=20)
#plt.title(r'$a = 0.5$')


plt.subplots_adjust(left = 0.07, right = 0.98, wspace = 0.17)
plt.savefig('ho2.pdf')



###############################
# end of heat maps


plt.show()




