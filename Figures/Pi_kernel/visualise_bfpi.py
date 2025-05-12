# Plot elements of Pi or the gradient of Pi in 3D
import numpy as np
import plotly.graph_objects as go

# Define a grid in 3D space centered at the origin
nx = 150
ny = 150
nz = 150

vx = 50
vy = 75
vz = 75
Xarray = np.linspace(-vx, vx, nx)
Yarray = np.linspace(-vy, vy, ny)
Zarray = np.linspace(-vz, vz, nz)


def functional(x, i,j, y, z, surval=0):
    r = np.sqrt(x**2 + y**2 + z**2)

    if   i ==1:
        v1 = x
    elif i ==2:
        v1 = y
    else:
        v1 = z
    if   j ==1:
        v2 = x
    elif j ==2:
        v2 = y
    else:
        v2 = z
    val = - 3*(v1*v2)  / r**5
    if i==j:
        val += 1/(r**3)
    return val - surval


def gradpi(i,j, k, x, y, z):
    r = np.sqrt(x**2 + y**2 + z**2)

    if   i ==1:
        v1 = x
    elif i ==2:
        v1 = y
    else:
        v1 = z

    if   j ==1:
        v2 = x
    elif j ==2:
        v2 = y
    else:
        v2 = z

    if   k ==1:
        v3 = x
    elif k ==2:
        v3 = y
    else:
        v3 = z

    val = 15*(v1*v2*v3) / r**7


    if   i==j:
        val -= 3 * v3 / r**5
    elif k==i:
        val -= 3 * v2 / r ** 5
    elif j==k:
        val -= 3 * v1 / r ** 5

    return val





threshold = 150e-7
coords = []
for xvals in Xarray:
    for yvals in Yarray:
        for zvals in Zarray:

            v = functional(xvals, i=3, j=3, y=yvals, z=zvals)
            if v > np.abs(threshold):
                coords.append([xvals, yvals, zvals, v])

#
# #THRESHOLD FOR GRADIENT OF PI
# threshold = 150e-6
# coords = []
# for xvals in Xarray:
#     for yvals in Yarray:
#         for zvals in Zarray:
#
#             v = gradpi(i=2, j=2, k=3, x=xvals, y=yvals, z=zvals)
#             if v > np.abs(threshold):
#                 coords.append([xvals, yvals, zvals, v])



pcloud = np.array(coords)

logpcloud = np.log10(pcloud[:, 3])


fig = go.Figure(data=go.Scatter3d(
    x=pcloud[:,0],
    y=pcloud[:,1],
    z=pcloud[:,2],
    mode='markers',
    marker=dict(size=3,
                color=np.log10(pcloud[:,3]),
                colorscale='Viridis',
                cmin=np.min(logpcloud),  # or choose a fixed range
                cmax=np.max(logpcloud)-2,
                )
))



fig.show()
