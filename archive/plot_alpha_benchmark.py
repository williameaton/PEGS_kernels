# The idea here will be to try and calculate the convolution of the gravity kernels directly
import ensightreader
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from funcs import EnData
from analytical import analytical, analytical_pot_homo_fullsphere
from matplotlib.gridspec import GridSpec
def norm(x):
    return x/np.max(np.abs(x))


R = 6371000
Rcmb = R*(1- 2891/6371)
plot_connectivity = False
G = 6.6723e-11

# Testing for intermediate kernel vector [1, 0, 0]
k1 = 1
k2 = 0
k3 = 0


nprocs = 100
fpath ="./NEX192/"

fig = plt.figure(layout="constrained", figsize=(10,4))
gs = GridSpec(8,9, figure=fig)

ax1 = fig.add_subplot(gs[:7,:3])
ax2 = fig.add_subplot(gs[:7,3:6])
ax3 = fig.add_subplot(gs[:7,6:])

# colorbars
ax1cb = fig.add_subplot(gs[-1,2:7])


X = np.array([])
Y = np.array([])
Z = np.array([])
V = np.array([])

for i_proc in range(nprocs):


    DProc = EnData(iproc=i_proc, fpath=fpath)
    DProc.load_ProcData('grav1_kernel')

    g1 = DProc.vardata["grav1_kernel"]

    mask = np.logical_and(DProc.node_coordinates[:,2] < 2, DProc.node_coordinates[:,2] > -2)
    x = DProc.node_coordinates[mask,0]
    y = DProc.node_coordinates[mask,1]
    z = DProc.node_coordinates[mask,2]
    v = g1[mask]
    X = np.append(X,x[:])
    Y = np.append(Y,y[:])
    Z = np.append(Z,z[:])
    V = np.append(V,v[:])

X*=R
Y*=R
Z*=R
from copy import copy
V0 = copy(V)


ff,aa = plt.subplots()
aa.scatter((X**2 + Y**2 + Z**2)**0.5, V * (6371000 * (5515.13*np.pi*6.67e-11)) )

plt.show()

analy  = analytical_pot_homo_fullsphere(dimrho=2600, R=R, x=X,y=Y,z=Z)

"""px = np.arange(1.5, 1.8, 0.001)
py = []
for p in px:
    diff = np.abs(V[:]/p - analy[:])
    py.append(np.sum(diff)/len(diff))

f, ax = plt.subplots()
ax.plot(px,py)
plt.show()"""

max = np.max(np.abs(V))

min = -max


diff     = V[:]-analy[:]
diffnorm     = V[:]/np.max(np.abs(V)) - analy[:]/np.max(np.abs(analy))


# Create colormap:
cmap = mpl.cm.RdBu_r
norm = mpl.colors.Normalize(vmin=min, vmax=max)

fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
             cax=ax1cb, orientation='horizontal', label='')


im1 = ax1.scatter(X, Y, c=analy, marker='o', cmap=cmap, vmin=min, vmax=max, s=5)
im2 = ax2.scatter(X, Y, c=V,     marker='o', cmap=cmap, vmin=min, vmax=max, s=5)
im3 = ax3.scatter(X, Y, c=diff,  marker='o', cmap=cmap, vmin=min, vmax=max, s=5)


for a in [ax1, ax2, ax3]:
    buffer = 1.05
    a.set_xlim([-R*buffer,R*buffer])
    a.set_ylim([-R*buffer,R*buffer])

    a.tick_params(
        axis='x',  # changes apply to the x-axis
        which='both',  # both major and minor ticks are affected
        bottom=False,  # ticks along the bottom edge are off
        top=False,  # ticks along the top edge are off
        labelbottom=False)

    a.axes.yaxis.set_ticklabels([])
    a.set_yticks([])

ax1.set_title('Analytical solution')
ax2.set_title('Real solution')
ax3.set_title('Difference')

# Create an array of theta values for the circle
theta = np.linspace(0, 2 * np.pi, 100)
x = Rcmb * np.cos(theta)
y = Rcmb * np.sin(theta)
for atmp in [ax1, ax2, ax3]:
    atmp.plot(x, y, 'k')

    atmp.spines[['right', 'top', 'bottom', 'left']].set_visible(False)

plt.show()
