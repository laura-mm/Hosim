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

#crit = np.loadtxt('crit_20_0.txt', delimiter=",")
#crih = np.loadtxt('crit_5_0.txt', delimiter=",")

g0 = np.loadtxt('measfp_2_0_0.txt', delimiter = ",")
g0 = g0.reshape(len(g0)/6, 6)
g0 = np.transpose(g0)
g1 = np.loadtxt('measfp_2_0_1.txt', delimiter = ",")
g1 = g1.reshape(len(g1)/6, 6)
g1 = np.transpose(g1)
g2 = np.loadtxt('measfp_2_0_2.txt', delimiter = ",")
g2 = g2.reshape(len(g2)/6, 6)
g2 = np.transpose(g2)

g = [g0, g1, g2]

meas = np.loadtxt('homeas_2_0.txt', delimiter=",")
meas = meas.reshape(3, grid + 1, 6)
meas = np.transpose(meas, (0, 2, 1)) # gam, meas, sig

col = ['m', 'c', 'y']
colo = ['m.', 'c.', 'y.']
lab = ['$\gamma = -\frac{1}{p-1}$', '$\gamma = 0$', '$\gamma = 1$']
mark = ['o', 's', '*', 'x', '+']
#lett = ['$(a)$', '$(b)$', '$(c)$', '$(d)$', '$(e)$', '$(f)$']
x = np.linspace(-1, 1, grid + 1)
s = np.power(10, x)

# sigma, z1, phi, M, q, help for g

# phi, M, q, diversity, dsq, h, for meas

# phi, M, diversity
###############################################################

plt.subplots(3, 1, figsize=(4,12))

plt.subplot(3, 1, 1) # phi
for gi in range(3):
	plt.semilogx(g[gi][0], g[gi][2], col[gi])
	plt.semilogx(s, meas[gi][0], colo[gi], marker = mark[gi], label = lab[gi])
	#plt.axvline(x = crit[gi], linestyle = '--', color = col[gi])
plt.xlim([10**(-1), 10**1])
#plt.ylabel(r'$\phi$')
plt.ylim([-0.2, 1.01])
#plt.text(6, 0.85, lett[0], fontsize=20)
#plt.title(r'$p = 2$')

plt.subplot(3, 1, 2) # M
for gi in range(3):
	plt.semilogx(g[gi][0], g[gi][3], col[gi])
	plt.semilogx(s, meas[gi][1], colo[gi], marker = mark[gi], label = lab[gi])
	#plt.axvline(x = crit[gi], linestyle = '--', color = col[gi])
plt.xlim([10**(-1), 10**1])
#plt.ylabel(r'$M^*$')
plt.ylim([0.2, 2])
#plt.text(6, 1.8, lett[2], fontsize=20)

plt.subplot(3, 1, 3) # diversity
for gi in range(3):
	plt.semilogx(g[gi][0], g[gi][3]*g[gi][3]/g[gi][4], col[gi])
	plt.semilogx(s, meas[gi][3], colo[gi], marker = mark[gi], label = lab[gi])
	#plt.axvline(x = crit[gi], linestyle = '--', color = col[gi])
plt.xlim([10**(-1), 10**1])
#plt.xlabel(r'$\sigma$')
#plt.ylabel(r'$\text{diversity }\frac{M^2}{q}$')
plt.ylim([0.3, 1])
#plt.text(6, 0.9, lett[4], fontsize=20)
#plt.legend(bbox_to_anchor=(0, -0.4, 2.17, 0.2), loc='lower left', ncol=5, mode="expand", borderaxespad=0)

#plt.subplots_adjust(left=0.09, right = 0.97, bottom = 0.13, top = 0.97, wspace = 0.17)
plt.savefig('hofive2.pdf')


plt.show()

