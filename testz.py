# fiveplots for higher order
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
grid = 10
col = ['m', 'c', 'y']
colo = ['m.', 'c.', 'y.']
lab = ['$\gamma = -\frac{1}{p-1}$', '$\gamma = 0$', '$\gamma = 1$']
mark = ['o', 's', '*', 'x', '+']
#lett = ['$(a)$', '$(b)$', '$(c)$', '$(d)$', '$(e)$', '$(f)$']
x = np.linspace(-1, 1, grid + 1)
s = np.power(10, x)

# sigma, z1, phi, M, q, help for measfp
# phi, M, q, diversity, dsq, h, for meas


pi = 1
#plt.subplots(5, 1)
mi = 1
gi = 2
g = np.loadtxt("Tmeasfp_{0}_{1}_{2}.txt".format((pi + 2), 20*(mi - 1), gi), delimiter = ",")
g = g.reshape(len(g)/6, 6)
g = np.transpose(g)
"""	
plt.subplot(5, 1, 1) # phi
plt.plot(g[1], g[2], col[gi])
#plt.xlim([-1, 1])
plt.ylabel("phi")
#plt.ylim([-0.2, 1.01])
"""
#plt.subplot(5, 1, 2) # M
plt.plot(g[1], g[0], col[gi])
plt.xlim([-1, 1])
plt.ylabel("M")
plt.ylim([-1, 1])
"""
plt.subplot(5, 1, 3) # sigma
plt.plot(g[1], g[0], col[gi])
#plt.xlim([-1, 1])
plt.ylabel("sigma")
#plt.ylim([-0.2, 1.01])

plt.subplot(5, 1, 4) # q
plt.plot(g[1], g[4], col[gi])
#plt.xlim([-1, 1])
plt.ylabel("q")
#plt.ylim([-0.2, 1.01])

plt.subplot(5, 1, 5) # help
plt.plot(g[1], g[5], col[gi])
#plt.xlim([-1, 1])
plt.ylabel("help")
#plt.ylim([-0.2, 1.01])
		
plt.savefig("testp3g2.pdf")

"""
plt.show()

