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
grid = 20
col = ['m', 'c', 'y']
colo = ['m.', 'c.', 'y.']
lab = ['$\gamma = -\frac{1}{p-1}$', '$\gamma = 0$', '$\gamma = 1$']
mark = ['o', 's', '*', 'x', '+']
#lett = ['$(a)$', '$(b)$', '$(c)$', '$(d)$', '$(e)$', '$(f)$']
x = np.linspace(-1, 1, grid + 1)
s = np.power(10, x)

# sigma, z1, phi, M, q, help for measfp
# phi, M, q, diversity, dsq, h, for meas

# phi, M, diversity
###############################################################
for pi in range (1):
	
	plt.subplots(3, 2, figsize=(12,8))

	for mi in range (2):

		meas = np.loadtxt("homeas_{0}_{1}.txt".format((pi + 2), 20*(mi - 1)), delimiter=",")
		meas = meas.reshape(3, grid + 1, 6)
		meas = np.transpose(meas, (0, 2, 1)) # gam, meas, sig
	
		crit = np.loadtxt("crit_{0}_{1}.txt".format((pi + 2), 20*(mi - 1)), delimiter=",")
			
		measfp = []
		for gi in range(3):
			g = np.loadtxt("measfp_{0}_{1}_{2}.txt".format((pi + 2), 20*(mi - 1), gi), delimiter = ",")
			g = g.reshape(len(g)/6, 6)
			g = np.transpose(g)
			measfp.append(g)
		
		plt.subplot(3, 2, 1 + mi) # phi
		for gi in range(3):
			plt.semilogx(measfp[gi][0], measfp[gi][2], col[gi])
			plt.semilogx(s, meas[gi][0], colo[gi], marker = mark[gi], label = lab[gi])
			plt.axvline(x = crit[gi], linestyle = '--', color = col[gi])
		plt.xlim([10**(-1), 10**1])
		if mi == 0:
			plt.ylabel("phi")
		#plt.ylabel(r'$\phi$')
		#plt.ylim([-0.2, 1.01])
		#plt.text(6, 0.85, lett[0], fontsize=20)
		plt.title("mu = " + str((2*mi) - 2))
		#plt.title(r'$p = 2$')
		
	
		plt.subplot(3, 2, 3 + mi) # M
		for gi in range(3):
			plt.semilogx(measfp[gi][0], measfp[gi][3], col[gi])
			plt.semilogx(s, meas[gi][1], colo[gi], marker = mark[gi], label = lab[gi])
			plt.axvline(x = crit[gi], linestyle = '--', color = col[gi])
		plt.xlim([10**(-1), 10**1])
		#plt.ylabel(r'$M^*$')
		if mi == 0:
			plt.ylabel("M")
		plt.ylim([0, 10])
		#plt.text(6, 1.8, lett[2], fontsize=20)
	
		plt.subplot(3, 2, 5 + mi) # diversity
		for gi in range(3):
			plt.semilogx(measfp[gi][0], measfp[gi][3]*measfp[gi][3]/measfp[gi][4], col[gi])
			plt.semilogx(s, meas[gi][3], colo[gi], marker = mark[gi], label = lab[gi])
			plt.axvline(x = crit[gi], linestyle = '--', color = col[gi])
		plt.xlim([10**(-1), 10**1])
		#plt.xlabel(r'$\sigma$')
		plt.xlabel("sigma")
		#plt.ylabel(r'$\text{diversity }\frac{M^2}{q}$')
		if mi == 0:
			plt.ylabel("diversity")
		#plt.ylim([0, 1])
		#plt.text(6, 0.9, lett[4], fontsize=20)
		#plt.legend(bbox_to_anchor=(0, -0.4, 2.17, 0.2), loc='lower left', ncol=5, mode="expand", borderaxespad=0)

	#plt.subplots_adjust(left=0.09, right = 0.97, bottom = 0.13, top = 0.97, wspace = 0.17)
	plt.savefig("hofive{0}.pdf".format(pi + 2))


plt.show()

