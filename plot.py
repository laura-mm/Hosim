# for plotting trajectories from simulations
import numpy as np
import matplotlib.pyplot as plt
N = 200
T = 2002
x = np.loadtxt('trajectoriesx.txt', delimiter=",")
x = x.reshape(N,len(x)/N)
y = np.loadtxt('trajectoriesy.txt', delimiter=",")
y = y.reshape((N,len(y)/N))

xF = np.loadtxt('trajectoriesxF.txt', delimiter=",")
xF = xF.reshape(N,len(xF)/N)
yF = np.loadtxt('trajectoriesyF.txt', delimiter=",")
yF = yF.reshape((N,len(yF)/N))

print(x.shape)

plt.figure(1)
for i in range(N):
	plt.plot(x[i])
#plt.tick_params(axis = x, pad = 2000000)
#plt.ylim([-1,10])
plt.figure(2)
for i in range(N):
	plt.plot(y[i])
plt.figure(3)
for i in range(N):
	plt.plot(yF[i])
plt.figure(4)
for i in range(N):
	plt.plot(yF[i])
plt.show()

