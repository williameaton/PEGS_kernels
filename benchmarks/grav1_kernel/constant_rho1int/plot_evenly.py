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

fig, ax = plt.subplots(figsize=(5,7))

wanted_sliced = []


if plot_connectivity:
    figmesh = plt.figure(figsize=(5, 7))
    axmesh = figmesh.add_subplot(111, projection='3d')



total = 0



X = np.array([])
Y = np.array([])
Z = np.array([])
V = np.array([])

for i_proc in range(nprocs):
    DProc = ProcData(proc=i_proc, fpath=fpath)


    DProc.load_var('grav1_kernel')
    g1 = DProc.data["grav1_kernel"]

    mask = np.logical_and(DProc.coord[:,2] < 0.000000000005, DProc.coord[:,2] > -0.000000000005)

    if  (mask==False).all():
        pass
    else:
        wanted_sliced.append(i_proc)


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
             orientation='horizontal', label='')


im1 = ax.scatter(X, Y, c=analy, marker='o', cmap=cmap, vmin=min, vmax=max, s=5)


for a in [ax]:
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


# Create an array of theta values for the circle
theta = np.linspace(0, 2 * np.pi, 100)
x = Rcmb * np.cos(theta)
y = Rcmb * np.sin(theta)
for atmp in [ax]:
    atmp.plot(x, y, 'k')

    atmp.spines[['right', 'top', 'bottom', 'left']].set_visible(False)



np.savetxt(fname='wanted_slices.txt', X=wanted_sliced)