import numpy as np
import matplotlib.pyplot as plt
from matplotlib import pyplot
from mpl_toolkits.mplot3d import Axes3D
from pylab import figure


# Load:
c = np.loadtxt('sem_node_coords', delimiter=',')

id = c[:,0]
x  = c[:,1]
y  = c[:,2]
z  = c[:,3]


fig = figure()
ax = fig.add_subplot(projection='3d')

ax.scatter(x, y, z, color='b')

for i in range(len(id)):
    ax.text(x[i], y[i], z[i], '%s' % (str(id[i])), size=5, zorder=1, color='k')

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')


plt.show()