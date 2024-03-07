import numpy as np
import matplotlib.pyplot as plt
from copy import copy

R = 1
Rcmb = (1-2891/6371)*R

# Create a mesh with spacing dx
nn = 61

x = np.linspace(-R, R, nn)
dx = x[1]-x[0]
dv = dx**3


# Create mesh:
X,Y,Z = np.meshgrid(x,x,x)

circle_mask = np.logical_and((X**2 + Y**2 + Z**2)**0.5 >= Rcmb, (X**2 + Y**2 + Z**2)**0.5 <= R)

x = X[circle_mask].flatten()
y = Y[circle_mask].flatten()
z = Z[circle_mask].flatten()

data = copy(x)*0

N = len(x)

# Define k
k11 = 1; k12=0; k13=0
k21 = 0; k22=0; k23=0
k31 = 0; k32=0; k33=1


count = 0
for master_elmt in range(N):
    sum = 0

    # Coords of master element:
    cx = x[master_elmt]
    cy = y[master_elmt]
    cz = z[master_elmt]

    if cz >= -0.00001 and cz <= 0.00001:

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
                r1 = cx - cxtmp
                r2 = cy - cytmp
                r3 = cz - cztmp

                mag = (r1**2 + r2**2 + r3**2)**0.5
                i_mag3 = 1/(mag**3)
                i_mag5 = 1/(mag**5)


                val =    k11*(i_mag3 - 3*(r1*r1)*i_mag5) +\
                         k22*(i_mag3 - 3*(r2*r2)*i_mag5) +\
                         k33*(i_mag3 - 3*(r3*r3)*i_mag5) +\
                         (k12+k21)*(- 3*(r1*r2)*i_mag5) +\
                         (k13+k31)*(- 3*(r1*r3)*i_mag5) +\
                         (k23+k32)*(- 3*(r2*r3)*i_mag5)









                #val =  6.6723e-11 * dv * (k1*r1 + k2*r2 + k3*r3) / ((r1**2 + r2**2 + r3**2)**1.5)
                sum += val * 6.6723e-11 * dv

        count += 1
        if count%50 ==0:
            print(f'Completed {count}')
        print(sum)
        data[master_elmt] = sum
    else:
        data[master_elmt] = np.nan



print(np.nanmin(data)); print(np.nanmax(data));



# Now plot:
fig = plt.figure()
ax = fig.add_subplot(111)#, projection='3d')
ax.scatter(x, y, c=data, cmap='RdBu_r')
ax.set_xlim([-1,1])
ax.set_ylim([-1,1])
plt.show()


np.savetxt(fname=f'fd_data_{nn}_{k11}_{k12}_{k13}_{k21}_{k22}_{k23}_{k31}_{k32}_{k33}', X=np.array([x,y,z,data]).T)

