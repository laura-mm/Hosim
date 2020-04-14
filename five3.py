# fiveplots for p = 3
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
s = np.linspace(0, 0.5, grid + 1)

# sigma, z1, phi, M, q, help for measfp
# phi, M, q, diversity, dsq, h, for meas

# phi, M, diversity
###############################################################

	
plt.subplots(3, 1)



meas = np.loadtxt("homeas_3_+15.txt", delimiter=",")
meas = meas.reshape(grid + 1, 6)
meas = np.transpose(meas) # meas, sig

crit = np.loadtxt("crit_3_15.txt", delimiter=",")

g = np.loadtxt("measfp_3_15.txt", delimiter = ",")
g = g.reshape(len(g)/6, 6)
g = np.transpose(g)


plt.subplot(3, 1, 1) # phi
plt.plot(g[0], g[2])
plt.scatter(s, meas[0], marker = '.')
plt.axvline(x = crit, linestyle = '--')
plt.axvline(x = 0.205543, linestyle = '--', color = 'k')

#plt.xlim([10**(-2), 10**1])

	
plt.subplot(3, 1,2) # M
plt.plot(g[0], g[3], color = 'r')
plt.scatter(s, meas[1], marker = '.')
plt.axvline(x = crit, linestyle = '--')
plt.axvline(x = 0.205543, linestyle = '--',  color = 'k')
#plt.xlim([10**(-2), 10**1])
#plt.ylim([0, 10])
	
plt.subplot(3, 1, 3) # diversity
plt.plot(g[0], g[3]*g[3]/g[4])
plt.scatter(s, meas[3], marker = '.')
plt.axvline(x = crit, linestyle = '--')
plt.axvline(x = 0.205543, linestyle = '--',  color = 'k')

#plt.xlim([10**(-2), 10**1])
#plt.ylim([0, 1])
plt.savefig("five3.pdf")

plt.figure(2)
plt.plot(g[1], g[2], marker = '.')

plt.figure(3)
plt.plot(g[1], g[0])

plt.figure(4)
plt.plot(g[1], g[3])
plt.ylim([-10, 10])

plt.figure(5)
plt.plot(g[1], g[4])
plt.ylim([-10, 10])

plt.show()

