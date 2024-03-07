import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
fig, ax = plt.subplots(3, 6, figsize=(14,7))

pimap = ['11', '22', '33', '12', '13', '23']

# Define receiver location:


# Create a grid:
n = 200
R = 6371
buffer = 1.05
x_r = [0,0,R]
xlin = np.linspace(-R*buffer, R*buffer, n)
ylin = np.linspace(-R*buffer, R*buffer, n)
zlin = np.linspace(-R*buffer, R*buffer, n)
Xd, Yd, Zd = np.meshgrid(xlin,ylin,zlin)

"""Xd = Xd.flatten()
Yd = Yd.flatten()
Zd = Zd.flatten()"""

# Create Radius mask:
rmask = (Xd**2 + Yd**2 + Zd**2)**0.5 > R

#Xd=Xd[rmask]
#Yd=Yd[rmask]
#Zd=Zd[rmask]

#maskval = 0.015
#xplotmask = np.abs(Xd) < maskval
#yplotmask = np.abs(Yd) < maskval
#zplotmask = np.abs(Zd) < maskval

X12 = x_r[0] - Xd
Y12 = x_r[1] - Yd
Z12 = x_r[2] - Zd
R12 = (X12**2 + Y12**2 + Z12**2)**0.5

# Calculate PI:
# pi = delta_ij / R3    -   3 (r_i - r_i')(r_j - r_j') / R5

pi = []
pi.append(1/(R12**3) - 3*(X12*X12)/(R12**5)) # 1, 11
pi.append(1/(R12**3) - 3*(Y12*Y12)/(R12**5)) # 2, 22
pi.append(1/(R12**3) - 3*(Z12*Z12)/(R12**5)) # 3, 33

pi.append(-3*(X12*Y12)/(R12**5)) # 4, 12
pi.append(-3*(X12*Z12)/(R12**5)) # 5, 13
pi.append(-3*(Y12*Z12)/(R12**5)) # 6, 23


for i in range(6):
    pi[i][rmask] = -999999

# Write to ENSIGHT to study:
filename='Pi.nc'
f = nc.Dataset(filename, 'w', format='NETCDF4')

# Create the dimensions:
x_dim = f.createDimension('x_dim', n)
y_dim = f.createDimension('y_dim', n)
z_dim = f.createDimension('z_dim', n)

# Creating the variables:
x   = f.createVariable('x', 'f4', ('x_dim',))
y   = f.createVariable('y', 'f4', ('y_dim',))
z   = f.createVariable('z', 'f4', ('z_dim',))

v_pi_11  = f.createVariable('pi_11', 'f4', ('x_dim', 'y_dim', 'z_dim',))
v_pi_22  = f.createVariable('pi_22', 'f4', ('x_dim', 'y_dim', 'z_dim',))
v_pi_33  = f.createVariable('pi_33', 'f4', ('x_dim', 'y_dim', 'z_dim',))
v_pi_12  = f.createVariable('pi_12', 'f4', ('x_dim', 'y_dim', 'z_dim',))
v_pi_13  = f.createVariable('pi_13', 'f4', ('x_dim', 'y_dim', 'z_dim',))
v_pi_23  = f.createVariable('pi_23', 'f4', ('x_dim', 'y_dim', 'z_dim',))

# Assigning values to the variables:
x[:] = xlin
y[:] = ylin
z[:] = zlin

v_pi_11[:,:,:] = pi[0]
v_pi_22[:,:,:] = pi[1]
v_pi_33[:,:,:] = pi[2]
v_pi_12[:,:,:] = pi[3]
v_pi_13[:,:,:] = pi[4]
v_pi_23[:,:,:] = pi[5]

print('Data written to file ', filename)
f.close()


# Calc max values of any pi:
mins = []
maxes = []
for i in range(6):
    mins.append(np.min(pi[i][pi[i] > -999998]))
    maxes.append(np.max(pi[i] ))

mins = np.array(mins)
maxes = np.array(maxes)

print(np.min(mins))
print(np.max(maxes))


"""
# Plotting:
for i in range(6):
    size = 0.7
    # XY plane
    ax[0,i].scatter(Xd[zplotmask], Yd[zplotmask], s=size, c=np.abs(pi[i][zplotmask]))

    # XZ plane
    ax[1,i].scatter(Xd[yplotmask], Zd[yplotmask], s=size, c=np.abs(pi[i][yplotmask]))

    # YZ plane
    ax[2,i].scatter(Yd[xplotmask], Zd[xplotmask], s=size, c=np.abs(pi[i][xplotmask]))


for row in range(3):
    for i in range(6):
        atmp = ax[row,i]

        #atmp.plot(xcirc,ycirc, 'k')
        atmp.set_xlim([-1.05, 1.05])
        atmp.set_ylim([-1.05, 1.05])

        atmp.tick_params(labelleft=False, labelbottom=False, left=False, bottom=False)
        atmp.spines['top'].set_visible(False)
        atmp.spines['right'].set_visible(False)
        atmp.spines['left'].set_visible(False)
        atmp.spines['bottom'].set_visible(False)

        if row==0:
            atmp.set_title(rf'$\Pi {str(pimap[i])}$')
            atmp.plot(x_r[0], x_r[1], 'r*')
        elif row==1:
            atmp.plot(x_r[0], x_r[2], 'r*')
        else:
            atmp.plot(x_r[1], x_r[2], 'r*')

ax[0,0].set_ylabel('XY plane')
ax[1,0].set_ylabel('XZ plane')
ax[2,0].set_ylabel('YZ plane')
"""


plt.show()
