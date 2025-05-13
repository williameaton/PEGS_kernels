import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append(f"../../../classes")
sys.path.append(f"../../../Figures")
from scipy import integrate
from finite_difference import FD, FdStation
from plot_all_stations import plot_all_stations
from wetools.plotting import save_figs_to_single_pdf, setup_we_mpl, subplots_hide_xaxes, align_y_labels
from colour_schemes import hex
setup_we_mpl()



# Load the strain data:
spfm_data_path = f'../../../strain/data/specfem/MDJ_kernel_forward_sim/'


fd = FD(from_readme=True, readme_fname='STRAIN_MDJ_README' )
a = fd.stations[0].snms = ['MDJ_33']


chls =  ['NN', 'EE', 'ZZ', 'ZN', 'EN', 'ZE']
fd.load_specfem_data(types=['strain'], dpath=spfm_data_path, samp_code='M')
ichl = 0
istn = fd.stations[0]

data = [ ]
for chl in chls:
    axi = 0
    d = istn.data
    time = d['time']

    # Integrate twice with respect to time

    chldata = d[istn.name + '_' + str(istn.cstn)][f'{chl[::-1]}.STRAIN']


    int1 = integrate.cumtrapz(y=chldata, x=time, initial=0)
    int2 = integrate.cumtrapz(y=int1, x=time,    initial=0)
    ichl +=1
    data.append(int2)


# Now we have the 6 channels that we need we want to compute the contraction H : H
fig, ax = plt.subplots(4)

figstrain, axstrain = plt.subplots()

ddata = np.array(data)
misfit = ddata[0,:]*0
for c in [0, 1, 2, 3, 4, 5]:

    if c <= 2:
        mult = 1
    else:
        mult = 2

    # Plot misfit
    misfit += mult* (ddata[c,:]**2)

    # Plot strain

    ax[0].plot(time, ddata[c,:])
    axstrain.plot(time, ddata[c,:])
    ax[0].set_ylabel('Strain')
    axstrain.set_ylabel('Strain')

    # Plot squared value
    ax[1].plot(time, ddata[c,:]**2)
    ax[1].set_ylabel('Squared h components')

misfit *= 0.5


chi = integrate.cumtrapz(y=misfit, x=time, initial=0)

# Plot misfit before time integration
ax[2].plot(time, misfit)
ax[2].set_ylabel('Before time int.')

# Plot final misfit
ax[3].plot(time, chi)
ax[3].set_ylabel('Final misfit')

plt.show()