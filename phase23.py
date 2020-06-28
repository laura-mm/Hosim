from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
import math
import mpld3

"""
crit = np.loadtxt("crit23.txt", delimiter = ",")
crit = crit.reshape(1001, 1001, 4) # mu2 mu3 (mu2, mu3, sig2, sig3)
crit = np.transpose(crit, (2, 0, 1)) # (mu2, mu3, sig2, sig3), mu2, mu3

mu2 = crit[0]
mu3 = crit[1]
sig2 = crit[2]
sig3 = crit[3]

#div = np.loadtxt("div23.txt", delimiter = ",")

div = div.reshape(1111, 6) # (mu2 & z1) (mu2, mu3, z1, sigtot, sig2, sig3)
div = np.transpose(div)

mu2d = div[0]
mu3d = div[1]
sig2d = div[4]
sig3d = div[5]

"""




divlines = [] # mu2, z1, (mu2, mu3, z1, sigtot, sig2, sig3)
fig1 = plt.figure(1)
ax1 = fig1.gca(projection = '3d')
#ax1.w_zaxis.set_scale('log')
fig2 = plt.figure(2)
ax2 = fig2.gca(projection = '3d')
#ax2.zaxis.set_scale('log')
for i in range(101):
	

	divline = np.loadtxt("div23_{0}.txt".format(i), delimiter = ",")
	divline = divline.reshape(len(divline)//6, 6) # z1, (mu2, mu3, z1, sigtot, sig2, sig3)
	divline = np.transpose(divline) # (mu2, mu3, z1, sigtot, sig2, sig3), z1
	divlines.append(divline)
	
	critline = np.loadtxt("crit23_{0}.txt".format(i), delimiter = ",")
	critline = critline.reshape(len(critline)//6, 6) # mu3, (mu2, mu3, sig2, sig3, maxsig, M)
	critline = np.transpose(critline) # (mu2, mu3, sig2, sig3, maxsig, M), mu3
	
	crit2 = np.loadtxt("critp2_{0}.txt".format(i), delimiter = ",")
	"""
	if (i == 10):
		print("mu2 = {0}".format(critline[0][0]))
		plt.figure(7)
		plt.plot([-4, crit2[1]], [crit2[2], crit2[2]], "r")
		plt.plot(divline[1], divline[4], "b")
		plt.figure(8)
		plt.plot(divline[1], divline[5], "b")
		plt.plot(critline[1], critline[4], "b")
		plt.plot(critline[1], critline[3], "r")
	"""
	#ax1.plot(critline[0], critline[1], critline[2], "r")
	ax1.plot([crit2[0], crit2[0]], [-4, crit2[1]], [crit2[2], crit2[2]], "r")
	ax1.plot(divline[0], divline[1], divline[4], "b")
	#ax1.plot(critline[0], critline[1], critline[5], "g")
	
	ax2.plot(divline[0], divline[1], divline[5], "b")
	ax2.plot(critline[0], critline[1], critline[4], "b")
	ax2.plot(critline[0], critline[1], critline[3], "r")
	#ax.plot(divline[0], divline[1], divline[4])
#ax.plot_trisurf(divline[0], divline[1], divline[5])
ax1.set_xlabel("mu2")
ax1.set_ylabel("mu3")
ax1.set_zlabel("sig2")
ax1.set_xlim(-4,6)
ax1.set_ylim(-4,6)
ax1.set_zlim(0,10)

ax2.set_xlabel("mu2")
ax2.set_ylabel("mu3")
ax2.set_zlabel("sig3")
ax2.set_xlim(-4,6)
ax2.set_ylim(-4,6)
ax2.set_zlim(0,10)



"""

fig3 = plt.figure(3)
ax3 = fig3.gca(projection = '3d')

fig4 = plt.figure(4)
ax4 = fig4.gca(projection = '3d')

fig5 = plt.figure(5)
ax5 = plt.axes()

fig6 = plt.figure(6)
ax6 = plt.axes()

i = 4

# mu2, mu3, (z1, mu2, mu3, sig2, M), z1


for j in range(11):
	s2t = np.loadtxt("s2test_{0}_{1}.txt".format(i,j), delimiter = ",")
	s2t = s2t.reshape(len(s2t)//5, 5) 
	s2t = np.transpose(s2t) 
	
	if (j == 0):
		print("mu2 = {0}".format(s2t[1][0]))
		
	if (j == 4):
		print("mu3 = {0}".format(s2t[2][0]))
		ax5.plot(s2t[0], s2t[3], "r")
		ax5.plot(s2t[0], s2t[4], "b")
		ax6.plot(s2t[3], s2t[4], "g")


	#ax3.plot(s2t[2], s2t[0], s2t[3], "r")
	ax3.plot(s2t[2], s2t[0], s2t[4], "b")
	ax4.plot(s2t[2], s2t[3], s2t[4], "g")

ax3.set_ylim(-50, 0)
ax3.set_zlim(-10, 10)
ax4.set_zlim(-10, 10)
#ax5.set_ylim(-10, 10)
#ax6.set_ylim(-10, 10)

"""

"""
test = np.loadtxt("s2test.txt", delimiter = ",")
test = test.reshape(len(test)//5, 5)
test = np.transpose(test)
print("mu2 = {0}".format(test[1][0]))
print("mu3 = {0}".format(test[2][0]))
plt.figure(1)
plt.plot(test[0], test[3])
plt.plot(test[0], test[4])
plt.ylim([-10,10])
plt.figure(2)
plt.plot(test[3], test[4])
plt.ylim([0,10])
plt.xlim([0,10])
"""

plt.show()
#mpld3.show()
