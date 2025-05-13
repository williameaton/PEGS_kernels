# Finite difference to compute STRAIN based on grav potential at numerous stations
# Note, only doing this for the surface gradients since the vertical gradients seem pretty good

# Requires phi measurements labelled by the format:
#               N
#              53
#              43
#W     31  32  33  34  35     E
#              23
#              13
#               S

# Data dir:
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import cumtrapz, trapz


def second_time_int(a, dx):
    a = cumtrapz(a, dx=dx, initial=0, axis=0)
    a = cumtrapz(a, dx=dx, initial=0, axis=0)
    return a

separation = 0.05
stn_lb = 'H'
first_n_timesteps = 10000000


ddir   = f'benchmark_homogenous_spheres/specfem/compare_FDs_step/data/{stn_lb}stns/raw'
outdir = f'benchmark_homogenous_spheres/specfem/compare_FDs_step/data/{stn_lb}stns/FD'

for stnno in range(1,7):
    stn = f'{stn_lb}S{stnno}'

    chlfreq = 'B'


    # Spatial separation - NOTE NOT INTERESTED IN THE FLATTENING HERE
    dEW = separation * (np.pi*6371000)/180 # 0.05 degrees of longitude
    dNS = separation * (np.pi*6371000)/180 # 0.05 degrees of latitude

    # Load phi data:
    phi = []

    # Load a 5 x 5 matrix of phi values for each station rows = SW and columns = WE
    # increasing index means more North or more East
    for i in range(1,6):
        phiEW = []
        for j in range(1,6):
            # West to East
            phiEW.append(np.loadtxt(f"{ddir}/{stn}_{i}{j}.YY.{chlfreq}XG.sem.ascii")[:first_n_timesteps,1])

        phi.append(phiEW)


    # Compute first derivatives along the stencil
    dp_dSN = []
    dp_dWE = [ ]
    for p in range(5):
        # del phi/ del y (NS)
        # at the EW index P
        dp_dSN.append( (-phi[4][p] + 8*phi[3][p]  - 8*phi[1][p] + phi[0][p])/(12*dNS)    )

        # del phi/ del x (EW)
        # at the NS index P
        dp_dWE.append( (-phi[p][4] + 8*phi[p][3] - 8*phi[p][1] + phi[p][0])/(12*dEW)    )


    # Compute second (cross derivatives) - they should be the same?
    #    -- del^2 phi/ del N del E
    d2p_dWEdNS = (-dp_dWE[4] + 8*dp_dWE[3] - 8*dp_dWE[1] + dp_dWE[0])/(12*dNS)
    #    -- del^2 phi/ del E del N
    d2p_dNSdWE = (-dp_dSN[4] + 8*dp_dSN[3] - 8*dp_dSN[1] + dp_dSN[0])/(12*dEW)


    #    -- del^2 phi/ del N del N
    NN = (-phi[4][2] + 16*phi[3][2] - 30*phi[2][2] + 16*phi[1][2] - phi[0][2])/(12*dNS*dNS)
    #    -- del^2 phi/ del E del E
    EE = (-phi[2][4] + 16*phi[2][3] - 30*phi[2][2] + 16*phi[2][1] - phi[2][0])/(12*dEW*dEW)


    # To compute the vertical derivatives we need to use the PGRAV data in the vertical direction
    pgrav_WE = []
    pgrav_SN = []
    for q in range(1,6):
        pgrav_WE.append(np.loadtxt(f"{ddir}/{stn}_3{q}.YY.{chlfreq}XZC.PGRAV.sem.ascii")[:first_n_timesteps, 1]) # Increasing q eastwards
        pgrav_SN.append(np.loadtxt(f"{ddir}/{stn}_{q}3.YY.{chlfreq}XZC.PGRAV.sem.ascii")[:first_n_timesteps, 1]) # Increasing q northwards

    #    -- del^2 phi/ del Z del E
    d2p_dZdWE = (-pgrav_WE[4] + 8*pgrav_WE[3] - 8*pgrav_WE[1] + pgrav_WE[0])/(12*dEW)
    #    -- del^2 phi/ del Z del N
    d2p_dZdSN = (-pgrav_SN[4] + 8*pgrav_SN[3] - 8*pgrav_SN[1] + pgrav_SN[0])/(12*dNS)


    # Write outputs to file:
    time = np.loadtxt(f"{ddir}/{stn}_33.YY.{chlfreq}XG.sem.ascii")[:first_n_timesteps,0]

    outchlnames = ['NN', 'EE', 'EZ', 'NZ', 'NE']
    outchldata  = [NN, EE, d2p_dZdWE, d2p_dZdSN, d2p_dWEdNS]

    for iout in range(5):
        outstr = f"{outdir}/{stn}_33.YY.{chlfreq}X{outchlnames[iout]}.STRAIN.sem.ascii"
        np.savetxt(fname=outstr, X=np.array([time, outchldata[iout]]).T)
        print(f'Written to {outstr}')
