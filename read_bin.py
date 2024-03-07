from scipy.io import FortranFile
import numpy as np


lat = [-70, -8]
lon = [-70, -65]


def asSpherical(x, y, z):
    #takes list xyz (single coord)
    r       =  np.sqrt(x*x + y*y + z*z)
    lat     =  90 - np.arccos(z/r)*180/ np.pi #to degrees
    long    =  np.arctan2(y,x)*180/ np.pi
    return np.array([r,lat,long]).T




Nsample = 2430000
proc = 4


strproc = str(proc)
if len(strproc)==1:
    strproc = '0'+strproc

print('Read coordinates')
coord = np.loadtxt(f'reg1_proc{proc}_WEcoords')[:Nsample, :]


sph = asSpherical(coord[:,0], coord[:,1], coord[:,2])
mask_lat = np.logical_and(sph[:,1] >= lat[0], sph[:,1] <= lat[1])
mask_lon = np.logical_and(sph[:,2] >= lon[0], sph[:,2] <= lon[1])
mask = np.logical_and(mask_lon, mask_lat)


print('Read fortran kernel file')

b = FortranFile(filename=f'DATABASES_MPI/proc0000{strproc}_reg1_rho_kernel.bin', mode='r',)
b = b.read_reals(dtype=float)[:Nsample][mask]

print(proc, max(b))


finlen = len(b)


arr = np.zeros((int(finlen), 4))
arr[:,:3] = coord[mask,:]
arr[:,3] = b
np.savetxt(f"proc_{proc}.csv", arr,
              delimiter = ",")

