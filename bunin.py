import numpy as np
import matplotlib.pyplot as plt

col = ['r', 'm', 'b', 'c', 'g']
#lab = ['$\gamma = -1$', '$\gamma = -0.5$', '$\gamma = 0$', '$\gamma = 0.5$', '$\gamma = 1$']
lab = ['g = -1', 'g = -0.5', 'g = 0', 'g = 0.5', 'g = 1']


for i in range (5):

	z = np.loadtxt("bunin2_{0}.txt".format(5*(i - 2)), delimiter = ",")
	z = z.reshape(len(z)/3, 3)
	z = np.transpose(z)


	plt.plot(z[0], np.minimum(z[1], z[2]), col[i])
	plt.plot(z[0], z[2], col[i], label = lab[i])
	#plt.axhline(y = 0, linestyle = '--', color = 'r')

plt.xlim([-5, 3])
plt.ylim([0, 5])
plt.xlabel('mu')
plt.ylabel('sigma')
#plt.legend(bbox_to_anchor=(0, -0.4, 2.17, 0.2), loc='lower left', ncol=5, mode="expand", borderaxespad=0)
plt.legend(bbox_to_anchor=(0.5, -0.1), loc = 'upper center', ncol = 5)
plt.subplots_adjust(bottom = 0.17, top = 0.93)


plt.savefig('bunin.pdf')

plt.figure(2)
z = np.loadtxt("bunin3_0.txt", delimiter = ",")
z = z.reshape(len(z)/3, 3)
z = np.transpose(z)

plt.plot(z[0], z[1])
plt.plot(z[0], z[2])
plt.axhline(y = 0, linestyle = '--', color = 'r')
plt.ylim([0, 10])
plt.xlim([-5, 3])

plt.show()
