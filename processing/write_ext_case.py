# Generate a case file outside of the SPECFEM setup
# TODO: Add transience
import numpy as np
import os

suffix = {'grav1_kernel' : 'grav1kernel',
          'grav2_kernel' : 'grav2kernel',
          'rho_not_prime': 'rhonotprimekernel',
          'rho_kernel'   : 'rhokernel',
          'alpha_kernel' : 'alphakernel',
          'beta_kernel'  : 'betakernel',
          'kappa_kernel' : 'kappakernel',
          'mu_kernel'    : 'mukernel',
          }

# SPECIFY THE VARIABLES YOU WANT IN THE CASE FILE HERE
var_list = ['grav1_kernel',
            'grav1_kernel',
            'inter_rho1',
            'inter_rho2',
            'inter_rho3',
            'inter_rho4',
            'inter_rho5',
            'alpha_kernel',
            'beta_kernel',
            'kappa_kernel',
            'mu_kernel',
            'rho_kernel',
            'rho_not_prime',
            ]

region = 1 # 1 = mantle, 2 = OC, 3 = IC
fpath ="../kernels/MDJ_kernel/cases"


#proc_list = range(nprocs)
proc_list = np.loadtxt("wanted_slices/wanted_slices_pegs_tohoku.txt").astype(int)

for i_proc in proc_list:
    casefname = f"{fpath}/reg1_proc{i_proc}.case"
    f   = open(casefname, "w")

    f.write("FORMAT\n")
    f.write("type:  ensight gold\n\n")
    f.write("GEOMETRY\n")
    f.write(f"model:    reg{int(region)}_proc{i_proc}.geo\n\n")

    f.write("VARIABLE\n")

    for v in var_list:

        try:
            suf = suffix[v]
        except:
            suf = v
        f.write(f'scalar per node:          1 {v} reg{region}_proc{i_proc}.{suf}\n\n')

    f.close()

