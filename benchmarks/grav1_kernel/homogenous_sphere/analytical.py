import matplotlib.pyplot as plt
import numpy as np
from matplotlib.gridspec import GridSpec

R = 1
G = 6.67e-11
pi = np.pi
k =  1
k1 = 1
k2 = 1
k3 = 1

r = np.linspace(0, 2*R, 100000)

inmask  = r <= R
outmask = r >= R

# For phi inside of the solid sphere:
phi_in = (2/3)*pi*G*k*(3*(R**2) - r[inmask]**2)
phi_out = G*(4/3)*pi*(R**3)*k  /r[outmask]

# g = - nabla phi
g_in = (4/3)*pi*G*k*r[inmask]
g_out =  G*(4/3)*pi*(R**3)/ (r[outmask]**2)



fig, ax = plt.subplots(4)

ax[0].plot(r[inmask],phi_in, 'k')
ax[0].plot(r[outmask],phi_out, 'k')

ax[1].plot(r[inmask],g_in, 'k')
ax[1].plot(r[outmask],g_out, 'k')

ax[0].axvline(R)
ax[1].axvline(R)


# Now let us plot the direction vectors
n = 45
x =  np.linspace(-R, R, n)
y =  np.linspace(-R, R, n)
z = 0 #np.linspace(-R, R, n)
X,Y,Z = np.meshgrid(x,y,z)

circle = ((X**2) + (Y**2) + (Z**2))**0.5 <= R
x = X[circle].flatten()
y = Y[circle].flatten()
z = Z[circle].flatten()



r = (x**2 + y**2 + z**2)**0.5
theta = np.arccos(z/r)
phi  = np.sign(y)*np.arccos(x/(x**2 + y**2)**0.5)



constant = - (4/3)*pi*G
xmag = constant*r*np.sin(theta)*np.cos(phi)
ymag = constant*r*np.sin(theta)*np.sin(phi)
zmag = constant*r*np.cos(theta)

fig, ax = plt.subplots(1,3)



# Select only the ones on y = 0
ax[0].quiver(x,y, k1*xmag, xmag*0)#, ymag, angles='xy', scale=0.07)
ax[1].quiver(x,y, xmag*0, k2*ymag)#, ymag, angles='xy', scale=0.07)


# Comparing the divergences:
fig = plt.figure(layout="constrained", figsize=(11,4.5))
gs = GridSpec(9,9, figure=fig)
ax1 = fig.add_subplot(gs[:-1,:3])
ax2 = fig.add_subplot(gs[:-1,3:6])
ax3 = fig.add_subplot(gs[:-1, 6:])

axcb = fig.add_subplot(gs[-1,:])



# Load the fd values:
fd = np.loadtxt('fd_data_55')
nanmask =  np.invert(np.isnan(fd[:,-1]))
fd = fd[nanmask, :]
x = fd[:,0]
y = fd[:,1]
z = fd[:,2]
v = fd[:,3]

# Now plotting the analytical divergence:
div = -(4/3)*pi*G * (k1*x + k2*y + k3*z)


bound = np.max(np.abs(div))

import matplotlib as mpl
cmap = mpl.cm.RdBu_r
norm = mpl.colors.Normalize(vmin=-bound, vmax=bound)

fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
             cax=axcb, orientation='horizontal', label='')


an = ax1.scatter(x,y, c=div, cmap=cmap, vmin=-bound, vmax=bound)
ax1.set_title('Analytical')

re = ax2.scatter(x,y, c=v, cmap=cmap, vmin=-bound, vmax=bound)
ax2.set_title('Convolution cartesian')

diff = ax3.scatter(x,y, c=v-div, cmap=cmap, vmin=-bound, vmax=bound)
ax3.set_title('Difference')
#diffnorm = ax[3].scatter(x,y, c= v/np.max(np.abs(v)) - div/np.max(np.abs(div)), cmap=cmap, vmin=-bound, vmax=bound )


plt.show()