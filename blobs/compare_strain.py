import sys
from copy import copy
import numpy as np
import matplotlib.pyplot as plt
sys.path.append(f"../classes")
from finite_difference import FD


blob_data_dir = '../data/data_blobs/'


chls = ['NN', 'EE', 'ZN', 'EN', 'ZE']


fd1 = FD(from_readme=True, readme_fname='strain_readme')
fd1.load_simple_spfmx_strain(dpath=f'{blob_data_dir}/test_1_midpoint_blob/raw')

fd3 = FD(from_readme=True, readme_fname='strain_readme')
fd3.load_simple_spfmx_strain(dpath=f'{blob_data_dir}/test_3_corrected_source_blob/raw')


fd4 = FD(from_readme=True, readme_fname='strain_readme')
fd4.load_simple_spfmx_strain(dpath=f'{blob_data_dir}/test_4_source_red_blob/raw')


fd5 = FD(from_readme=True, readme_fname='strain_readme')
fd5.load_simple_spfmx_strain(dpath=f'{blob_data_dir}/test_5_central_cylinder/forward_altervpvsrho/raw')



fd6 = FD(from_readme=True, readme_fname='strain_readme')
fd6.load_simple_spfmx_strain(dpath=f"{blob_data_dir}/test_6_source_cylinder/raw")


fdPREM1C = FD(from_readme=True, readme_fname='strain_readme')
fdPREM1C.load_simple_spfmx_strain(dpath="../forward_simulations/Tohoku/data/PREM_1C/raw")



fdPREM2C = FD(from_readme=True, readme_fname='strain_readme')
fdPREM2C.load_simple_spfmx_strain(dpath="../forward_simulations/Tohoku/data/PREM_2C/raw")


# Load the strain directly for directory:

fig, ax = plt.subplots(fd1.nstns, 5, figsize=(16, 9.5))

ls = ['-', '-', '--', ':']

ichl = 0
for chl in chls:
    isim = 0
    for fff in [fdPREM1C, fdPREM2C, fd3]:
        axi = 0

        for istn in fff.stations:
            d = istn.data
            time = d['time']

            ldirect, = ax[axi, ichl].plot(time,  1e15 * d[f'{chl[::-1]}.STRAIN'], linestyle=ls[isim])

            ax[axi, 0].set_ylabel(istn.name)
            axi += 1
        isim += 1
    ichl += 1

plt.show()