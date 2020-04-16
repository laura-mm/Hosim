# for plotting trajectories from simulations
import numpy as np
import matplotlib.pyplot as plt
N = 20
x = np.loadtxt('trajectories.txt', delimiter=",")
x = x.reshape(N,len(x)/(3*N), 3)
x = np.transpose(x, (0, 2, 1))

#print(x.shape)

plt.figure(1)
for i in range(N):
	plt.plot(x[i][0], x[i][1])
#plt.tick_params(axis = x, pad = 2000000)
#plt.ylim([-1,10])
plt.figure(2)
for i in range(N):
	plt.plot(x[i][0], x[i][2])
plt.show()

