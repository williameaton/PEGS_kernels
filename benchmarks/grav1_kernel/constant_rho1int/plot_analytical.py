import matplotlib.pyplot as plt
import numpy as np


N = 500

G = 6.67e-11
pi = np.pi
Rcmb = 2890/6371
R    = 1.0


k1 = 1
k2 = 1
k3 = 1

X = np.linspace(-1, 1, N)
Y = np.linspace(-1, 1, N)
Z = np.linspace(-1, 1, N)

fig, ax = plt.subplots(figsize=(10.5,6.5))
fig.set_tight_layout(True)

z = 0

x,y = np.meshgrid(X,Y)
#x, y, z = np.meshgrid(X,Y,Z)
r = (x**2 + y**2 + z**2)**0.5


mask = np.logical_and(r>=Rcmb, r<=R)







c = ax.scatter(x[mask], y[mask], c=sol[mask], s =0.1, cmap=plt.cm.RdBu_r)

plt.colorbar(c)



for rad in [Rcmb,R]:
    theta = np.linspace(0, 2 * np.pi, 1000)
    x = rad * np.cos(theta)
    y = rad * np.sin(theta)
    ax.plot(x, y, 'k', linewidth=0.8)


ax.set_aspect('equal', adjustable='box')

ax.set_xlim([-1,1])
ax.set_ylim([-1,1])

ax.text(x=0.85, y=0.95, s=f'k1: {k1}')
ax.text(x=0.85, y=0.90, s=f'k2: {k2}')
ax.text(x=0.85, y=0.85, s=f'k3: {k3}')
ax.text(x=0.65, y=-0.95, s=f'Plane Z = {z}')

ax.set_xlabel('X Axis')
ax.set_ylabel('Y Axis')

plt.show()