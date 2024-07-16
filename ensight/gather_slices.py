# Script uses wanted slices text file to move necessary data to different directory
import numpy as np
import os

master_dir = f'./MDJ_KERNEL'

# Load wanted slices from text file
ws = np.loadtxt(f"{master_dir}/wanted_slices_MDJ.txt").astype(int)


#os.mkdir(f'{master_dir}/wanted_slices_dir')

ftypes = ['case', 'geo', 'inter_rho1', 'inter_rho2', 'inter_rho3', 'inter_rho4', 'inter_rho5', 'rhokernel', 'grav1kernel', 'grav2kernel']

for iproc in ws:
    for ftyp in ftypes:
        os.rename(f"{master_dir}/ensight_data/reg1_proc{iproc}.{ftyp}", f"{master_dir}/wanted_slices_dir/reg1_proc{iproc}.{ftyp}")
