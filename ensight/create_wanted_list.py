# The idea here will be to try and calculate the convolution of the gravity kernels directly
import ensightreader
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from load_proc_data import ProcData
from matplotlib.gridspec import GridSpec
import math as m


def cart2sph(x,y,z):
    XsqPlusYsq = x**2 + y**2
    r = (XsqPlusYsq + z**2)**0.5              # r
    elev = np.arctan2(z,XsqPlusYsq**0.5)     # theta
    az = np.arctan2(y,x)                           # phi
    return r, 180*elev/np.pi,  180*az/np.pi


R = 6371000
Rcmb = R*(1- 2891/6371)
G = 6.6723e-11


nprocs = 1536
fpath ="./MDJ_KERNEL/ensight_data"



lat_lim = np.array([25, 45])
lon_lim = np.array([105, 150])






wanted_sliced = []


X = np.array([])
Y = np.array([])
Z = np.array([])
V = np.array([])

for i_proc in range(nprocs):
    DProc = ProcData(proc=i_proc, fpath=fpath)

    DProc.load_var('grav1_kernel')

    xc = np.array(DProc.coord[:,0])
    yc = np.array(DProc.coord[:,1])
    zc = np.array(DProc.coord[:,2])

    r, lat, lon = cart2sph(xc, yc, zc) # - 180 to + 180, 90 to -90
    r*=R

    # Get and select coords on the surface
    rc = (xc**2 + yc**2 + zc**2)**0.5

    # Filter by radius close to surface
    surface = r==R
    latsurf = lat[surface]
    lonsurf = lon[surface]

    lat_mask = np.logical_and(latsurf > lat_lim[0], latsurf < lat_lim[1])
    lon_mask = np.logical_and(lonsurf > lon_lim[0], lonsurf < lon_lim[1])

    if (lon_mask==True).any() and (lat_mask==True).any():
        wanted_sliced.append(i_proc)
    else:
        pass


np.savetxt(fname='wanted_slices_Tohoku.txt', X=wanted_sliced)