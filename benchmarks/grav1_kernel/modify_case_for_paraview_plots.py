# The idea here will be to try and calculate the convolution of the gravity kernels directly
import numpy as np
import os

nprocs = 1536
fpath ="./plot_1st_kernel/"

for i_proc in range(nprocs):
    casefname = f"{fpath}/reg1_proc{i_proc}.case"


    file1   = open(casefname, "a")  # append mode


    file1.write(f'scalar per node:          1 grav1kernel reg1_proc{i_proc}.grav1kernel\n\n')

    file1.close()

