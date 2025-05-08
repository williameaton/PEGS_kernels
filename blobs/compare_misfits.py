# Compare this simulation against the 'normal' 1D prem model
import matplotlib.pyplot as plt
import sys
import numpy as np
sys.path.append('../classes')
from pegs import create_pegs
import pickle

blob_data_dir = '../data/data_blobs/'

forward_dir = '../data/data_forward_pegs'

# Load case without blobs
spec2c = create_pegs(f'{forward_dir}/PREM_2C/',
                   code='SPECFEM', three_component=True)
spec2c.load()

spec1c = create_pegs(f'{forward_dir}/PREM_1C/',
                   code='SPECFEM', three_component=True)
spec1c.load()

# Load case with blobs:
test1 = create_pegs(f'{blob_data_dir}/test_1_midpoint_blob/',           code='SPECFEM', three_component=True)
test2 = create_pegs(f'{blob_data_dir}/test_2_source_blob/'  ,           code='SPECFEM', three_component=True)
test3 = create_pegs(f'{blob_data_dir}/test_3_corrected_source_blob/',   code='SPECFEM', three_component=True)
test4 = create_pegs(f'{blob_data_dir}/test_4_source_red_blob/'  ,       code='SPECFEM', three_component=True)

test5 = create_pegs(f'{blob_data_dir}/test_5_central_cylinder/forward_altervpvsrho',       code='SPECFEM', three_component=True)
test6 = create_pegs(f'{blob_data_dir}/test_6_source_cylinder/'  ,       code='SPECFEM', three_component=True)



test1.plot_colour    = 'red'
test2.plot_colour   = 'green'



for c in [test1, test2, test3, test4, test5, test6]:
    c.load()
    c.cut_to_plot()

codes = [spec1c, test5, test6]

# Some plotting stuff:
figsize  = (6, 7.5)         # Figure size
lw_syn = 1                  # Synthetic data line width
grid_lines_on = False       # Grid lines?
vlinmin = 0.015             # Min y value of vertical lines
network = {'MAJO': 'IU', 'FUK': 'BO', 'SHR': 'BO', 'INU': 'G', 'BJT': 'IC', 'MDJ': 'IC', 'XAN': 'IC', 'INCN': 'IU',
           'MA2': 'IU', 'ULN': 'IU', 'NE93': 'YP'}

# Load the epicentral distances in km for annotation:
with open('../Figures/data_for_plots/epi_dist.pkl', 'rb') as fp:
    epidists = pickle.load(fp)


# Create figures and set tight layouts:
fig, ax     = plt.subplots(3, figsize=figsize)     # Figure for 1D Prem
fig.set_tight_layout(True)


# Stores line objects for mpl legend in loop
leg_lines = []
leg_lines_3D = []

# Loop over each station
ictr = 0
stn = 'MDJ'

# Plot each code for this station
codectr = 0

timemask = 164


fig, axmisfit = plt.subplots()

misfits = []

for icode in codes:
    mask = icode.time[stn] < timemask
    time = icode.time[stn][mask]

    misfit = icode.z_pegs[stn][mask] * 0
    iax = 0

    for chl in 'zne':

        chl_trace = icode.traces[f"{chl}_pegs"][stn][mask] * 1e-9

        ltmp, = ax[iax].plot(time, chl_trace, linewidth=lw_syn, label=icode.plot_label)

        ax[iax].set_ylabel(f'{chl.upper()} ' + "Pegs signal\n" + r"[$m$ $s^{-2}$]")

        # Add to misfit:
        misfit += chl_trace**2

        iax += 1

    axmisfit.plot(time, misfit, linewidth=lw_syn, label=icode.plot_label)


    # Misfit is half * integral
    # Now intgrate misfit over time to get misfit func:
    misint = 0.5 * np.trapz(x=time, y=misfit)
    if codectr == 0:
        print(f'PREM 1C: {misint }')
    elif codectr ==1:
        print(f'PREM 2C: {misint}')
    else:
        print(f'Test {codectr-1}: {misint}')

    misfits.append(misint)
    codectr += 1



axmisfit.set_xlabel('Time [s]')
axmisfit.set_ylabel('Misfit function $\chi$\n[$m^2$ $s^{-3}$]')

#axmisfit.set_xlim([120, 160])

#axmisfit.legend(['PREM 1C', 'PREM 2C', 'Test 1', 'Test 2', 'Test 3', 'Test 4'], ncols=2)

lbs = ['PREM 1C', '450 km midpoint', '75 km source']
axmisfit.legend(lbs, ncols=1)
ax[0].legend(lbs, ncols=1)
ax[-1].set_xlabel('Time [s]')


figsc, axsc = plt.subplots()

misfits = np.array(misfits)[[0, -2, -1]]

misfits = 100*(misfits - misfits[0])/misfits[0]
print(misfits)
axsc.scatter(np.arange(len(misfits)), np.array(misfits))
axsc.set_ylabel('% Variation from PREM 1C')

axsc.set_xticks(np.arange(len(misfits)))
axsc.set_xticklabels(lbs)

plt.show()



