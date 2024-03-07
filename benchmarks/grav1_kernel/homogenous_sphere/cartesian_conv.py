import numpy as np
import matplotlib.pyplot as plt
from copy import copy

# Define k
k1 = 1; k2=1; k3=1
R  = 1
Rcmb = (1-2891/6371)*R

# Create a mesh with spacing dx
nn = 55

x = np.linspace(-R, R, nn)
dx = x[1]-x[0]
dv = dx**3

# Create mesh:
X,Y,Z = np.meshgrid(x,x,x)

circle_mask =   (X**2 + Y**2 + Z**2)**0.5 <= R
x = X[circle_mask].flatten()
y = Y[circle_mask].flatten()
z = Z[circle_mask].flatten()

r = (x**2 + y**2 + z**2)**0.5
data = copy(r)*0
N = len(r)



zthreshold = 0.00001

count = 0
for master_elmt in range(N):
    sum = 0

    # Coords of master element:
    cx = x[master_elmt]
    cy = y[master_elmt]
    cz = z[master_elmt]

    if cz >= -zthreshold and cz <= zthreshold:

        # calculate the convolution:
        for i_elmt in range(N):
            # Ignore elements that are not in the circle, or are not in the mantle
            if i_elmt==master_elmt:
                pass
            else:
                # Local element coords:
                cxtmp = x[i_elmt]
                cytmp = y[i_elmt]
                cztmp = z[i_elmt]

                # caculate r' - r vector:
                r1 = cxtmp - cx
                r2 = cytmp - cy
                r3 = cztmp - cz



                val = dv * (k1*r1 + k2*r2 + k3*r3) / ((r1**2 + r2**2 + r3**2)**1.5)
                sum += val

        count += 1
        if count%50 ==0:
            print(f'Completed {count}')
        data[master_elmt] = sum*6.67e-11
    else:
        data[master_elmt] = np.nan


print(np.max(sum))
print(np.min(sum))

# Now plot:
fig = plt.figure()
ax = fig.add_subplot(111)#, projection='3d')
ax.scatter(x, y, c=data, cmap='RdBu_r')
ax.set_xlim([-1,1])
ax.set_ylim([-1,1])
plt.show()


np.savetxt(fname=f'fd_data_{nn}', X=np.array([x,y,z,data]).T)
