# Script uses wanted slices text file to move necessary data to different directory
import numpy as np
import os

master_dir = f'../blobs/midpoint_blob/ensight_setup'

orig_dir = f'density'
new_dir  = f'wanted'

# Load wanted slices from text file
#ws = np.loadtxt(f"{master_dir}/wanted_slices_MDJ.txt").astype(int)
ws = np.loadtxt(f"wanted_MDJ.txt").astype(int)


# #os.mkdir(f'{master_dir}/wanted_slices_dir')

#ftypes = ['case', 'geo', 'inter_rho1', 'inter_rho2', 'inter_rho3', 'inter_rho4', 'inter_rho5', 'rhokernel', 'grav1kernel', 'grav2kernel']
ftypes = [ 'density']

for iproc in ws:
    for ftyp in ftypes:

        deadname = f"{master_dir}/{orig_dir}/reg1_proc{iproc}.{ftyp}"
        newname  = f"{master_dir}/{new_dir}/reg1_proc{iproc}.{ftyp}"
        os.rename(deadname, newname)
        #print(f"{deadname} --> {newname}")
print('Done')