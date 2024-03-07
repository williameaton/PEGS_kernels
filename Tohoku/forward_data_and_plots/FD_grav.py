import numpy as np
import matplotlib.pyplot as plt

networks = {'MAJO':'IU', 'MDJ': 'IC','INU': 'GG',  'INCN': 'IU', 'NE93':'YP', 'BJT':'IC', 'MA2':'IU', 'XAN':'IC', 'ULN':'IU' }
stations = ['MAJO', 'INU', 'MDJ', 'INCN', 'NE93', 'BJT', 'MA2', 'XAN', 'ULN']


dir   = 'elastic_test/elastic_copy/OUTPUT_FILES'
FDdir = 'elastic_test/elastic_copy/FD_GRAV'

for stn in stations:
    net = networks[stn]
    separation = 2000

    phis = []
    # Load the potentials at each depth:
    for depth in [0, 2, 6, 8]:

        tmp = np.loadtxt(f"{dir}/{stn}_{depth}.{net}.MXG.sem.ascii")
        t   = tmp[:,0]
        phi = tmp[:,1]

        # Store the phis for each depth
        phis.append(phi)

    # Calculate the FD:

    grad = (-phis[0] + 8*phis[1] - 8*phis[2] + phis[3])/(12*separation)


    # Load the grav output:
    X = np.array([t, grad]).T
    np.savetxt(f"{FDdir}/{stn}_4.{net}.MXZC.PGRAV.sem.ascii", X=X)
