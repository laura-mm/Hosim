from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
import math
import plotly.express as px
import pandas as pd





#df = px.data.gapminder().query("country=='Brazil'")
#fig = px.line_3d(df, x="gdpPercap", y="pop", z="year")


"""
df = pd.DataFrame(dict(
    X=[0,1,2,3, 1,2,3,4], 
    Y=[0,2,3,1, 1,3,4,2], 
    Z=[0,3,1,2, 1,4,2,3],
    color=["a", "a", "a", "a", "b", "b", "b", "b"]
))



fig = px.line_3d(df, x='X', y='Y', z='Z', color="color")
fig.show()

"""

divlines = [] # mu2, z1, (mu2, mu3, z1, sigtot, sig2, sig3)
#fig1 = plt.figure(1)
#ax1 = fig1.gca(projection = '3d')
#fig2 = plt.figure(2)
#ax2 = fig2.gca(projection = '3d')
for i in range(1):
	

	divline = np.loadtxt("div23_{0}.txt".format(i), delimiter = ",")
	divline = divline.reshape(len(divline)//6, 6) # z1, (mu2, mu3, z1, sigtot, sig2, sig3)
	divline = np.transpose(divline) # (mu2, mu3, z1, sigtot, sig2, sig3), z1
	divlines.append(divline)
	
	critline = np.loadtxt("crit23_{0}.txt".format(i), delimiter = ",")
	critline = critline.reshape(len(critline)//6, 6) # mu3, (mu2, mu3, sig2, sig3, maxsig, M)
	critline = np.transpose(critline) # (mu2, mu3, sig2, sig3, maxsig, M), mu3
	
	crit2 = np.loadtxt("critp2_{0}.txt".format(i), delimiter = ",")
	"""
	df = pd.DataFrame(dict(
		X = divline[0], 
		Y = divline[1], 
		Z = divline[4]
	))
	fig = px.line_3d(df, x='X', y='Y', z='Z')
	"""
	
	fig = px.line_3d(divline[0], divline[1], divline[4])
	
	#ax1.plot([crit2[0], crit2[0]], [-4, crit2[1]], [crit2[2], crit2[2]], "r")
	#ax1.plot(divline[0], divline[1], divline[4], "b")
	
	#ax2.plot(divline[0], divline[1], divline[5], "b")
	#ax2.plot(critline[0], critline[1], critline[4], "b")
	#ax2.plot(critline[0], critline[1], critline[3], "r")
"""
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





fig.show()


