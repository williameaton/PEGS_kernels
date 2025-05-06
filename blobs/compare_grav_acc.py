# Compare this simulation against the 'normal' 1D prem model
import matplotlib.pyplot as plt
import sys
import numpy as np
sys.path.append('../classes')
from pegs import create_pegs
import pickle

# Load case without blobs
spec = create_pegs('../forward_simulations/Tohoku/data/globe_isotropic_prem/assuming_perfect_sphere/NEX256', code='SPECFEM')

blob_data_dir = '../data/data_blobs/'
# Load case with blobs:
test1 = create_pegs(f'{blob_data_dir}/test_1_midpoint_blob', code='SPECFEM', three_component=True)
test2 = create_pegs(f'{blob_data_dir}/test_2_source_blob', code='SPECFEM', three_component=True)
test3 = create_pegs(f'{blob_data_dir}/test_3_corrected_source_blob', code='SPECFEM', three_component=True)


test1.plot_colour   = 'red'
test2.plot_colour   = 'green'
test3.plot_colour   = 'purple'


codes = [test1, test2, test3]

for c in codes:
    c.load()
    c.cut_to_plot()

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
fig, ax     = plt.subplots(figsize=figsize)     # Figure for 1D Prem


for f in [fig]:
    f.set_tight_layout(True)

# Add grid lines
if grid_lines_on:
    for j in range(4):
        for a in [ax]:
            a.axvline(j*100, ymin=vlinmin, ymax=1,  color='k',linestyle='--', alpha=0.05)

# Stores line objects for mpl legend in loop
leg_lines = []
leg_lines_3D = []

# Loop over each station
ictr = 0
for stn in spec.stns:

    # Plot each code for this station
    for icode in codes:
        ltmp, = ax.plot(icode.time[stn], 9 + icode.z_grav[stn] - ictr, color=icode.plot_colour, linewidth=lw_syn, label=icode.plot_label)
        # Only add to legend string once (not for every station)
        if ictr == 0:
            leg_lines.append(ltmp)


    # Add station name and horizontal lines
    for a in [ax]:
        a.text(x=-130, y= 9 - ictr, s=stn, weight="bold", horizontalalignment='center', fontsize=10)
        a.text(x=-130, y= 9 - ictr-0.2, s=f"{int(epidists[stn])} km", weight="bold", horizontalalignment='center', fontsize=7.5)
        a.axhline(y=9 - ictr, xmin=0/450, xmax=100/450, color='grey', linewidth=1, alpha=0.3)



    ictr += 1



# Format legend
leg_titles    = ['SPECFEM', 'TEST 1', 'TEST 2', 'TEST 3']

legend   = ax.legend(leg_lines,leg_titles, fontsize=8)
for l in [legend]:
    frame = l.get_frame()
    frame.set_facecolor('white')
    frame.set_edgecolor('black')
    l.get_frame().set_alpha(None)

# Formatting the limits, labels etc
for a in [ax]:
    a.set_ylim([0, 9.5])
    a.set_xlim([-100, 350])
    a.axvline(0, ymin=vlinmin, ymax=1, linewidth=2, color='k')
    a.set_xlabel('Time since rupture initiation [s]', weight='bold')

    # Remove standard labels
    for xlabel_i in a.axes.get_yticklabels():
        xlabel_i.set_fontsize(10.0)
        xlabel_i.set_visible(False)

    # Scale bar for plot
    a.errorbar(x=300, y=2.5, xerr=0, yerr=0.5, color='k', linewidth=2, elinewidth=2, capsize=10 )
    a.text(x=305, y=2.4, s=r'1 $nm/s^2$' )

    # Remove ugly spines!
    a.tick_params(labelleft=False, left=False)
    a.spines['top'].set_visible(False)
    a.spines['right'].set_visible(False)
    a.spines['left'].set_visible(False)

    # X tick labels
    a.set_xticks([-100, 0, 100, 200, 300])
    a.set_xticklabels(('-100', '0', '100', '200', '300'))

plt.show()



