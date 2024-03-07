import numpy as np
from gll import gll


ngllx = 5
nglly = 5
ngllz = 5

# Get positions of GLL and weights
x, xw = gll(ngllx-1)
y, yw = gll(nglly-1)
z, zw = gll(ngllz-1)

x = np.array(x)
y = np.array(y)
z = np.array(z)



Lderiv = np.zeros((ngllx, ngllx))

# sampling L_j' at points x_p:
for p in range(ngllx):
    for j in range(ngllx):
        for l in range(ngllx):
            if j!=l:
                k = (1/(x[j]-x[l]))
                for m in range(ngllx):
                    if np.logical_and(m!=j, m!=l):
                        k *= (x[p] - x[m])/(x[j]-x[m])

                Lderiv[p,j] += k

