# Visualising the Pi(r, r') kernel
# In the adjoint source we require Pi_ij \hat{n}_j  where n is the z direction of the receiver coordinate system
# Let us therefore plot this in a receiver coordinated system, i.e. XYZ where Z is the up direction of the receiver,
# X and Y are East and North
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from wetools.plotting import setup_we_mpl
setup_we_mpl()

def norm(x):
    return x/np.max(np.abs(x))

# Create a sphere of radius 1:
theta, phi = np.mgrid[0:2 * np.pi:100j, 0.0001:0.1:100000j]
x = np.cos(theta) * np.sin(phi)
y = np.sin(theta) * np.sin(phi)
z = np.cos(phi)

# In our case the receiver is at 0,0,1
r =( (x**2) + (y**2) + (z-1)**2)**0.5


r3 = r**3
r5 = r**5


pi_n_1 = -3*(x-0)*(z-1)/r5
pi_n_2 = -3*(y-0)*(z-1)/r5
pi_n_3 = 1/r3
"""
# For colormap of the log, we want to round to the nearest whole? number
log3 = np.log10((pi_n_3))
vmin3 = np.floor(np.min(log3))
vmax3 = np.ceil(np.max(log3))
vspan = vmax3 - vmin3
log3 = (log3 - vmin3)/vspan

vrange = np.arange(vmin3, vmax3+1)
nints = len(vrange)

# Create custom colorbar
n = nints-1  # Number of discrete colors you want
colors = plt.cm.YlGnBu_r(np.linspace(0, 1, n))  # Sample n colors from a continuous colormap
discrete_cmap = LinearSegmentedColormap.from_list("custom_cmap", colors, N=n)

figz = plt.figure()
axz = figz.add_subplot( projection='3d')
zplot = axz.plot_surface(x, y, z, facecolors=discrete_cmap(log3))


mappable = plt.cm.ScalarMappable(cmap=discrete_cmap)
mappable.set_array(log3)
mappable.set_clim(0, 1)  # Normalized range [0, 1]

figz.colorbar(mappable, ax=axz)

cbax_z = figz.get_axes()[1]
cbax_z.set_yticks((vrange-vmin3)/vspan, vrange)"""


# Question is what is the closest GLL point, which determines the max amplitude.
# for the MDJ source it is found in an element
# We can see that the minimum distance in cartesian is given in variable 'dist'
gll_coords = np.loadtxt('./data_for_plots/Pi_kernel/MDJ_rec_coords', skiprows=1)
rec_coord = gll_coords[0,:]
dist = 100
for i in range(125):
    gll_dist = np.sum((gll_coords[i+1, :] - rec_coord)**2)**0.5
    if gll_dist < dist:
        dist = gll_dist

max_ir3 = 1/dist**3





# Load the number of GLL per radius:
ngll = np.loadtxt('./data_for_plots/Pi_kernel/MDJ_ngll_within_rad', skiprows=1)[:,1:]
ngll = ngll[ngll[:, 0].argsort()]



# Another way to think about the z plot is:
zc = np.linspace(1, 0.9, 10000000)
xc = (1 - zc**2)**0.5
yc = 0

rc =( (xc**2) + (zc-1)**2)**0.5

piz = 1/rc**3

fig, ax = plt.subplots()



# Sort the outputs from 1536 processors saying how many GLL are within threshold range
gll_dist = []
for i in range(1,101):
    thresh = i*0.004

    npw = np.where(np.abs(ngll[:,0]-thresh) < 0.00001)

    gll_dist.append(np.sum(ngll[npw,1]))
gll_dist = np.array(gll_dist)


ax.plot(rc, piz)
ax.plot(np.linspace(0.004, 0.404, 100), gll_dist)
ax.scatter(dist, max_ir3)


thresh  = 10**4
xthresh = thresh**(-1/3)
ax.axvline(xthresh, ymin=0, ymax=(np.log10(thresh)-1)/10 )
ax.axhline(thresh, xmin=0, xmax=xthresh/0.4)


ax.set_yscale('log')
ax.set_yscale('log')

ax.set_xlabel('Radial distance from source')
ax.set_xlim([0, 0.4])
ax.set_ylim([1e1, 1e11])
ax.legend([r'Magnitude decay of $1/r^3$', 'Number of GLL points on surface'])
plt.show()