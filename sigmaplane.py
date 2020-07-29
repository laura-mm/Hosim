from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
import math
import mpld3

#mu2val = [-4, 1, 6]
#mu3val = [-4, 0, 6]

mu2val = [-4.0, -2.0, 2.0, -4.0, -2.0, 2.0, -3.0, -2.0, 2.0]
mu3val = [5.0, 2.0, 2.0, 2.0, 0.0, 0.0, 1.0, -2.0, -2.0]

plt.subplots(3, 3)

for i1 in range(9):
	mu2 = mu2val[i1]
	mu3 = mu3val[i1]
		
	crit = np.loadtxt('planecrit_{0}.txt'.format(i1), delimiter = ',')
	crit = crit.reshape(len(crit)//4, 4)
	crit = np.transpose(crit) # z1, sig2, sig3, M
		
	smax = np.loadtxt('planemax_{0}.txt'.format(i1), delimiter = ',')
	smax = smax.reshape(len(smax)//5, 5)
	smax = np.transpose(smax) # z1, sig2, sig3, sigtot, M
		
	div = np.loadtxt('planediv_{0}.txt'.format(i1), delimiter = ',')
	div = div.reshape(len(div)//2, 2)
	div = np.transpose(div) # sig2, sig3
		
	plt.subplot(3, 3, i1 + 1)
	plt.plot(crit[1], crit[2])
	plt.plot(smax[1], smax[2])
	plt.plot(div[0], div[1])
	plt.xlim([0, 2])
	plt.ylim([0, 2])
		
		

plt.show()
