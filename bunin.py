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

plt.savefig('bunin2.pdf')

plt.figure(2)
z0 = np.loadtxt("bunin3_-5.txt", delimiter = ",")
z1 = np.loadtxt("bunin3_-4.txt", delimiter = ",")
z2 = np.loadtxt("bunin3_0.txt", delimiter = ",")
z3 = np.loadtxt("bunin3_5.txt", delimiter = ",")
z4 = np.loadtxt("bunin3_10.txt", delimiter = ",")
z = [z0, z1, z2, z3, z4]
m0 = np.loadtxt("mu3_-5.txt", delimiter = ",")
m1 = np.loadtxt("mu3_-4.txt", delimiter = ",")
m2 = np.loadtxt("mu3_0.txt", delimiter = ",")
m3 = np.loadtxt("mu3_5.txt", delimiter = ",")
m4 = np.loadtxt("mu3_10.txt", delimiter = ",")
m = [m0, m1, m2, m3, m4]
for i in range(5):

	z[i] = z[i].reshape(len(z[i])/3, 3)
	z[i] = np.transpose(z[i])
	plt.semilogy(z[i][0], z[i][2], col[i])
	m[i] = m[i].reshape(len(m[i])/3, 3)
	m[i] = np.transpose(m[i])
	plt.semilogy(m[i][0], m[i][1], col[i])
	#plt.plot(m[i][0], m[i][2], col[i])
	
plt.axhline(y = 0, linestyle = '--', color = 'k')
plt.ylim([-1, 2])
plt.xlim([-1, 1])





plt.savefig('signegg.pdf')

plt.show()
