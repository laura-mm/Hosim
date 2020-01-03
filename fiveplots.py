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

g0 = np.loadtxt('measfp_2_-20_0.txt', delimiter = ",")
g0 = g0.reshape(len(g0)/6, 6)
g0 = np.transpose(g0)
g1 = np.loadtxt('measfp_2_-20_1.txt', delimiter = ",")
g1 = g1.reshape(len(g1)/6, 6)
g1 = np.transpose(g1)
g2 = np.loadtxt('measfp_2_-20_2.txt', delimiter = ",")
g2 = g2.reshape(len(g2)/6, 6)
g2 = np.transpose(g2)

m0 = [g0, g1, g2] # mu = -2

g0 = np.loadtxt('measfp_2_0_0.txt', delimiter = ",")
g0 = g0.reshape(len(g0)/6, 6)
g0 = np.transpose(g0)
g1 = np.loadtxt('measfp_2_0_1.txt', delimiter = ",")
g1 = g1.reshape(len(g1)/6, 6)
g1 = np.transpose(g1)
g2 = np.loadtxt('measfp_2_0_2.txt', delimiter = ",")
g2 = g2.reshape(len(g2)/6, 6)
g2 = np.transpose(g2)

m1 = [g0, g1, g2] # mu = 0

g0 = np.loadtxt('measfp_2_20_0.txt', delimiter = ",")
g0 = g0.reshape(len(g0)/6, 6)
g0 = np.transpose(g0)
g1 = np.loadtxt('measfp_2_20_1.txt', delimiter = ",")
g1 = g1.reshape(len(g1)/6, 6)
g1 = np.transpose(g1)
g2 = np.loadtxt('measfp_2_20_2.txt', delimiter = ",")
g2 = g2.reshape(len(g2)/6, 6)
g2 = np.transpose(g2)

m2 = [g0, g1, g2] # mu = 2

m = [m0, m1, m2]



meas0 = np.loadtxt('homeas_2_-20.txt', delimiter=",") # mu = -2
meas0 = meas0.reshape(3, grid + 1, 6)
meas0 = np.transpose(meas0, (0, 2, 1)) # gam, meas, sig

meas1 = np.loadtxt('homeas_2_0.txt', delimiter=",") # mu = 0
meas1 = meas1.reshape(3, grid + 1, 6)
meas1 = np.transpose(meas1, (0, 2, 1)) # gam, meas, sig

meas2 = np.loadtxt('homeas_2_20.txt', delimiter=",") # mu = 2
meas2 = meas2.reshape(3, grid + 1, 6)
meas2 = np.transpose(meas2, (0, 2, 1)) # gam, meas, sig

meas = [meas0, meas1, meas2]

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

plt.subplots(3, 3, figsize=(12,12))

for mi in range (3):
	plt.subplot(3, 3, 1 + mi) # phi
	for gi in range(3):
		plt.semilogx(m[mi][gi][0], m[mi][gi][2], col[gi])
		plt.semilogx(s, meas[mi][gi][0], colo[gi], marker = mark[gi], label = lab[gi])
		#plt.axvline(x = crit[gi], linestyle = '--', color = col[gi])
	plt.xlim([10**(-1), 10**1])
	if mi == 0:
		plt.ylabel("phi")
	#plt.ylabel(r'$\phi$')
	plt.ylim([-0.2, 1.01])
	#plt.text(6, 0.85, lett[0], fontsize=20)
	plt.title("mu = " + str((2*mi) - 2))
	#plt.title(r'$p = 2$')
	

	plt.subplot(3, 3, 4 + mi) # M
	for gi in range(3):
		plt.semilogx(m[mi][gi][0], m[mi][gi][3], col[gi])
		plt.semilogx(s, meas[mi][gi][1], colo[gi], marker = mark[gi], label = lab[gi])
		#plt.axvline(x = crit[gi], linestyle = '--', color = col[gi])
	plt.xlim([10**(-1), 10**1])
	#plt.ylabel(r'$M^*$')
	if mi == 0:
		plt.ylabel("M")
	plt.ylim([0.2, 2])
	#plt.text(6, 1.8, lett[2], fontsize=20)

	plt.subplot(3, 3, 7 + mi) # diversity
	for gi in range(3):
		plt.semilogx(m[mi][gi][0], m[mi][gi][3]*m[mi][gi][3]/m[mi][gi][4], col[gi])
		plt.semilogx(s, meas[mi][gi][3], colo[gi], marker = mark[gi], label = lab[gi])
		#plt.axvline(x = crit[gi], linestyle = '--', color = col[gi])
	plt.xlim([10**(-1), 10**1])
	#plt.xlabel(r'$\sigma$')
	plt.xlabel("sigma")
	#plt.ylabel(r'$\text{diversity }\frac{M^2}{q}$')
	if mi == 0:
		plt.ylabel("diversity")
	plt.ylim([0.3, 1])
	#plt.text(6, 0.9, lett[4], fontsize=20)
	#plt.legend(bbox_to_anchor=(0, -0.4, 2.17, 0.2), loc='lower left', ncol=5, mode="expand", borderaxespad=0)

#plt.subplots_adjust(left=0.09, right = 0.97, bottom = 0.13, top = 0.97, wspace = 0.17)
plt.savefig('hofive2.pdf')


plt.show()

