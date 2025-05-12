# Edits all case files based on a single template case file:
import numpy as np
import os
# Define processor case files to edit
#procs = np.loadtxt("wanted_slices.txt").astype(int)
procs = np.loadtxt("wanted_MDJ.txt").astype(int)
#procs = np.arange(96)

# Input/output directories
#fpath     = "/Users/eaton/Documents/Princeton/PEGS_kernels/kernels/MDJ_kernel/ensight/"
#fpath_out = "/Users/eaton/Documents/Princeton/PEGS_kernels/kernels/MDJ_kernel/ensight"

fpath     = "/Users/eaton/Documents/Princeton/PEGS_kernels/kernels/MDJ_pegs_kernel/z_only_kernel"
fpath_out = "/Users/eaton/Documents/Princeton/PEGS_kernels/kernels/MDJ_pegs_kernel/z_only_kernel"

# Example processor:
ex_proc = 206

for p in procs:
    file1 = open(f'{fpath}/reg1_proc{ex_proc}.case', 'r')
    lines = file1.readlines()
    file1.close()
    file1 = open(f'{fpath_out}/reg1_proc{p}.case', 'w')
    # Replace line processor number:
    for i in range(len(lines)):
        if lines[i].find(f'proc{ex_proc}.')!=-1:
            lines[i] = lines[i].replace(f'proc{ex_proc}.', f'proc{p}.')
    file1.writelines(lines)
    file1.close()



