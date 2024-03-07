# The idea here will be to try and calculate the convolution of the gravity kernels directly
import ensightreader
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from load_proc_data import ProcData
from analytical import analytical, analytical_alpha
from matplotlib.gridspec import GridSpec

R = 6371000
Rcmb = R*(1- 2891/6371)
plot_connectivity = False
G = 6.6723e-11

k1 = 1
k2 = 1
k3 = 1

def norm(x):
    return x/np.max(np.abs(x))

write_values_to_case = True

nprocs = 864
fpath ="./NEX192/"

fig = plt.figure(layout="constrained", figsize=(10,7))
gs = GridSpec(15,9, figure=fig)

ax1 = fig.add_subplot(gs[:7,:3])
ax2 = fig.add_subplot(gs[:7,3:6])
ax3 = fig.add_subplot(gs[:7,6:])
# FD plots:
ax4 = fig.add_subplot(gs[7:14,:3])
ax5 = fig.add_subplot(gs[7:14,3:6])
ax6 = fig.add_subplot(gs[7:14,6:])

# colorbars
ax1cb = fig.add_subplot(gs[-1,2:7])


if plot_connectivity:
    figmesh = plt.figure(figsize=(5, 7))
    axmesh = figmesh.add_subplot(111, projection='3d')



total = 0



X = np.array([])
Y = np.array([])
Z = np.array([])
V = np.array([])

for i_proc in range(nprocs):
    print(i_proc)
    DProc = ProcData(proc=i_proc, fpath=fpath)


    DProc.load_var('grav1_kernel')
    g1 = DProc.data["grav1_kernel"]

    mask = np.logical_and(DProc.coord[:,2] < 0.000000000005, DProc.coord[:,2] > -0.000000000005)
    x = DProc.coord[mask,0]
    y = DProc.coord[mask,1]
    z = DProc.coord[mask,2]
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

V*= -1


analy = analytical(k1,k2,k3, Rcmb, X,Y,Z)


max = np.max(np.abs(analy))

min = -max
diff     = V[:]-analy[:]



# Create colormap:
cmap = mpl.cm.RdBu_r
norm = mpl.colors.Normalize(vmin=min, vmax=max)

fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
             cax=ax1cb, orientation='horizontal', label='')


im1 = ax1.scatter(X, Y, c=analy, marker='o', cmap=cmap, vmin=min, vmax=max, s=5)
im2 = ax2.scatter(X, Y, c=V,     marker='o', cmap=cmap, vmin=min, vmax=max, s=5)
im3 = ax3.scatter(X, Y, c=diff,  marker='o', cmap=cmap, vmin=min, vmax=max, s=5)


# Plot the FD:
fd = np.loadtxt(f'fd_data_61_{k1}_{k2}_{k3}')
# get non-nan values:
nanmask =  np.invert(np.isnan(fd[:,-1]))
fd = fd[nanmask, :]

x = fd[:,0]
y = fd[:,1]
z = fd[:,2]
fdv = fd[:,-1]
analyFD = analytical(k1,k2,k3, Rcmb,x,y,z)


diff = analyFD - fd[:,-1]

im4 = ax4.scatter(x, y, c=fdv,   marker='o', cmap=cmap, s=5, vmin=min, vmax=max)
im5 = ax5.scatter(x, y, c=analyFD,         marker='o', cmap=cmap, s=5, vmin=min, vmax=max)
im6 = ax6.scatter(x, y, c=diff,      marker='o', cmap=cmap, s=5, vmin=min, vmax=max)


for a in [ax1, ax2, ax3, ax4, ax5, ax6]:
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
ax4.set_title('Convolution')
ax5.set_title('Analytical at Grid points')
ax6.set_title('Difference')


# Create an array of theta values for the circle
theta = np.linspace(0, 2 * np.pi, 100)
x = Rcmb * np.cos(theta)
y = Rcmb * np.sin(theta)
for atmp in [ax1, ax2, ax3, ax4, ax5, ax6]:
    atmp.plot(x, y, 'k')

    atmp.spines[['right', 'top', 'bottom', 'left']].set_visible(False)

plt.show()
#plt.savefig('grav1_benchmark.pdf', format='pdf')

"""
if plot_connectivity:
    valpha = 0.3
    i=0
    axmesh.plot([cx[i+0], cx[i+1], cx[i+2], cx[i+3], cx[i+0]],
                [cy[i+0], cy[i+1], cy[i+2], cy[i+3], cy[i+0]],
                [cz[i+0], cz[i+1], cz[i+2], cz[i+3], cz[i+0]], 'k', alpha=valpha)
    i=4
    axmesh.plot([cx[i+0], cx[i+1], cx[i+2], cx[i+3], cx[i+0]],
                [cy[i+0], cy[i+1], cy[i+2], cy[i+3], cy[i+0]],
                [cz[i+0], cz[i+1], cz[i+2], cz[i+3], cz[i+0]], 'k', alpha=valpha)
    for i in range(4):
        axmesh.plot([cx[i+0], cx[i+4]],
                    [cy[i+0], cy[i+4]],
                    [cz[i+0], cz[i+4]], 'k', alpha=valpha)"""