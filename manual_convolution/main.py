# The idea here will be to try and calculate the convolution of the gravity kernels directly
import ensightreader
import matplotlib.pyplot as plt
import numpy as np
from gll import gll
from load_proc_data import ProcData
plot_connectivity = False

G = 6.67e-11

ngllx = 5
nglly = 5
ngllz = 5

fpath ="/Users/eaton/Downloads/download/download/"

fig = plt.figure(figsize=(5, 7))
ax = fig.add_subplot(111)

if plot_connectivity:
    figmesh = plt.figure(figsize=(5, 7))
    axmesh = figmesh.add_subplot(111, projection='3d')

nprocs = 864
major_proc = 12


total = 0

# First we need to load in the major proc stuff - this will be the one we are actually calculating the values for:
MainProc =  ProcData(proc=major_proc, fpath=fpath)




X = np.array([])
Y = np.array([])
Z = np.array([])
V = np.array([])

for i_proc in range(nprocs):
    if i_proc!= major_proc:
        pass
    else:

        DProc = ProcData(proc=i_proc, fpath=fpath)
        # load the variable at top of page
        DProc.load_var("rho1_int1")
        DProc.load_var("rho1_int2")
        DProc.load_var("rho1_int3")

        # For whole processor:
        kp1 = DProc.data["rho1_int1"]
        kp2 = DProc.data["rho1_int2"]
        kp3 = DProc.data["rho1_int3"]

        print('NEED TO CHANGE JACOBIAN TO REAL VALUE')
        # Load jacobians
        DProc.load_var('rho_kernel')



        # Loop over each element and calculate its kernel value
        for ielmt in range(DProc.part.number_of_elements):



            # Connectivity of single element
            conn = DProc.connectivity[ielmt] -1 # python starts at 0!
            # Coorinates of single element
            cx = DProc.coord[conn,0]
            cy = DProc.coord[conn,1]
            cz = DProc.coord[conn,2]

            eval1 = kp1[conn]
            eval2 = kp2[conn]
            eval3 = kp3[conn]

            rminrd =

            # only corner elements are available
            for i in range(8):
                dotprod = eval1(i)*

                total +=  (1/3)**3 * dotprod jacobian * lagrange_poly






total *= G


"""
        mask = np.logical_and(coord[:,1] < 0.02, np.logical_and(coord[:,1] > -0.02, coord[:,0]>0))
        x = DProc.coord[mask,0]
        y = DProc.coord[mask,1]
        z = DProc.coord[mask,2]
        v = d[mask]

        X = np.append(X,x[:])
        Y = np.append(Y,y[:])
        Z = np.append(Z,z[:])
        V = np.append(V,v[:])


vrange = 1e-15

im = ax.scatter( X, Z, c=V, marker='o', vmin=-vrange, vmax=vrange, cmap=plt.cm.RdBu_r, s=5)

ax.set_xlim([0,1])
ax.set_ylim([-1,1])
#ax.set_zlim([-1,1])

fig.colorbar(im, ax=ax)


plt.show()
"""



"""
            if plot_connectivity:
                valpha = 0.3
                i=0
                axmesh.plot([cx[i+0], cx[i+1], cx[i+2], cx[i+3], cx[i+0]],
                            [cy[i+0], cy[i+1], cy[i+2], cy[i+3], cy[i+0]],
                            [cz[i+0], cz[i+1], cz[i+2], cz[i+3], cz[i+0]], 'k', alpha=valpha)
                i=4
                axmesh.plot([cx[i+0], cx[i+1], cx[i+2], cx[i+3], cx[i+0]],
                            [cy[i+0], cy[i+1], cy[i+2], cy[i+3], cy[i+0]],
                            [cz[i+0], cz[i+1], cz[i+2], cz[i+3], cz[i+0]], 'k', alpha=valpha)
                for i in range(4):
                    axmesh.plot([cx[i+0], cx[i+4]],
                                [cy[i+0], cy[i+4]],
                                [cz[i+0], cz[i+4]], 'k', alpha=valpha)"""