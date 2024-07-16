# Edits all case files based on a single template case file:
import numpy as np

# Define processor case files to edit
procs = np.arange(96) #np.loadtxt("wanted_slices.txt").astype(int)

# Input/output directories
fpath     = "./STRAIN/adding_surface_src_to_specfem/benchmark_src/fordown/"
fpath_out = f"./STRAIN/adding_surface_src_to_specfem/benchmark_src/fordown/"

# Example processor:
ex_proc = 0

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
