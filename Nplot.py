import numpy as np
import matplotlib.pyplot as plt
import math


muval = []
muval.append([-1.0, 0.0])
muval.append([-0.25, 0.0])
gammaval = []
gammaval.append([-1.0, 0.0, 1.0])
gammaval.append([-0.5, 0.0, 1.0])
start = [-1.0, -1.5]
# N, phi, M, q, diversity, dsq, h, otherphi
data = np.loadtxt("Nplot.txt", delimiter = ",")
data = data.reshape(2, 2, 3, 21, 20, 8)
data = np.transpose(data, (0, 1, 2, 3, 5, 4))



mint3 = [0, 2]
gint3 = [0, 3, 5]

crit = [] # gi, mi
for mi in range(2):
	c = np.loadtxt("crit3_{0}.txt".format(mint3[mi]), delimiter = ",")
	crit.append(c)
crit = np.transpose(crit)




for pi in range(2):
	#pi = pi + 1
	p = pi + 2
	
	plt.subplots(2, 3, figsize = (18, 12))
	figno = 0
	plt.suptitle("p = {0}".format(p))
	
	for mi in range(2):
		#mi = mi + 1
		mu = muval[pi][mi]
		
		if (pi == 0):
			crit2 = np.loadtxt("crit_{0}_{1}.txt".format((pi + 2), 10*(mi - 1)), delimiter=",")
		for gi in range (3):
			figno = figno + 1
			#gi = gi + 2
			gamma = gammaval[pi][gi]
			M = []
			S = []
			for si in range(21):
				#plt.figure(figno)
				sigma = 10**(0.1*si + start[pi])
				S.append(sigma)
				x = np.log(data[pi][mi][gi][si][0])
				y = np.log(data[pi][mi][gi][si][0]*data[pi][mi][gi][si][1])
				#plt.plot(x, y, "x")
				m, b = np.polyfit(x, y, 1)
				M.append(m)
				#plt.plot(x, m*x + b)
				#plt.title("p = {0}, mu = {1}, gamma = {2}, sigma = {3}".format(p, mu, gamma, sigma))
			plt.subplot(2, 3, figno)
			plt.semilogx(S, M, "x")
			if (pi == 0):
				plt.axvline(x = crit2[gi], linestyle = '--')
			if (pi == 1):
				plt.axvline(x = crit[gint3[gi]][mi], linestyle = '--')
			if (gi == 0):
				plt.ylabel("mu = {}".format(mu))
			if (mi == 0):
				plt.title("gamma = {0}".format(gamma))
			if (mi == 1):
				plt.xlabel("sigma")
	plt.savefig("Nplot{}.pdf".format(p))
plt.show()
