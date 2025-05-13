# Plots Tohoku PEGS from SPECFEM compared with QSSP, AXITRA and the real data
import matplotlib.pyplot as plt
import sys
import numpy as np
sys.path.append('../../classes')
from pegs import create_pegs
import pickle

plot_grav = True

# Some plotting stuff:
leg_titles = ['MINEOS', 'QSSP', 'AXITRA', 'SPECFEM']
figsize  = (10, 7.5)         # Figure size
lw_syn = 2                   # Synthetic data line width
vlinmin = 0.015              # Min y value of vertical lines




# Load the epicentral distances in km for annotation:
with open('../../Figures/data_for_plots/epi_dist.pkl', 'rb') as fp:
    epidists = pickle.load(fp)



# Load specfem for perfect sphere
spec = create_pegs('data/globe_isotropic_prem/assuming_perfect_sphere/NEX256', code='SPECFEM')
spec.load()

# Load QSSP
qssp = create_pegs('data/qssp/original_AK135/output/', code='QSSP')
qssp.load()
qssp.cut_to_plot()

# Load AXITRA
axitra = create_pegs('data/Vallee_data/synthetic/', code='AXITRA', scale=1)
axitra.load()
axitra.cut_to_plot()

# Create figures and set tight layouts:
fig, axes     = plt.subplots(1,2, figsize=figsize)  # Figure for 1D Prem
fig.set_tight_layout(True)

# Add grid lines
if grid_lines_on:
    for j in range(4):
        for i in range(2):
            axes[i].axvline(j*100, ymin=vlinmin, ymax=1,  color='k', linestyle='--', alpha=0.05)

# Stores line objects for mpl legend in loop
leg_lines    = []

axitra.plot_colour = 'purple'

# Loop over the grav and zacc:
for iax in range(2):
    if iax == 0:
        plot_grav = True
        axtitle = 'Gravitational acceleration'
    else:
        plot_grav = False
        axtitle = 'Elastic acceleration'
    ax = axes[iax]

    # Loop over each station
    ictr = 0
    for stn in spec.stns:


        # Decided to add the Normal Modes Juhel code data:
        # This code was provided by Kevin Juhel
        # set PEGS and time vector
        network = {'MAJO': 'IU', 'FUK': 'BO', 'SHR': 'BO', 'INU': 'G', 'BJT': 'IC', 'MDJ': 'IC', 'XAN': 'IC', 'INCN': 'IU',
                   'MA2': 'IU', 'ULN': 'IU', 'NE93': 'YP'}
        nm_zacc = np.load('data/Juhel_2019_NM/KJ19/PEGS.' + network[stn] + '.' + stn + '.Z.proc_acc.npy')
        nm_pegs = np.load('data/Juhel_2019_NM/KJ19/PEGS.' + network[stn] + '.' + stn + '.Z.pegs_proc.npy')
        if plot_grav:
            plotvar = nm_pegs - nm_zacc
        else:
            plotvar = - nm_zacc
        time     = np.arange(0, nm_pegs.size, 1.0) - 100
        timemask = np.logical_and(time < spec.cutoff[stn], time >= 0)
        lnm, = ax.plot(time[timemask], 9 - ictr + (10 ** 9 * plotvar[timemask]), color='#BB5566', linewidth=lw_syn)
        if ictr == 0:
            leg_lines.append(lnm)






        # Plot each code for this station
        for icode in [qssp, axitra, spec]:

            if plot_grav:
                plotvar = icode.z_grav[stn]
            else:
                plotvar = icode.z_acc[stn]
            ltmp, = ax.plot(icode.time[stn], 9 + plotvar - ictr, color=icode.plot_colour, linewidth=lw_syn)
            # Only add to legend string once (not for every station)
            if ictr == 0:
                leg_lines.append(ltmp)

        # For left-most axis:
        if iax == 0:
            # Add station name and horizontal lines
            ax.text(x=-50, y= 9 - ictr, s=stn, weight="bold", horizontalalignment='center', fontsize=10)
            ax.text(x=-50, y= 9 - ictr-0.2, s=f"{int(epidists[stn])} km", weight="bold", horizontalalignment='center', fontsize=7.5)


        ictr += 1


    ax.set_title(axtitle, weight='bold')

    # Format legend
    legend   = ax.legend(leg_lines, leg_titles, fontsize=8)
    for l in [legend]:
        frame = l.get_frame()
        frame.set_facecolor('white')
        frame.set_edgecolor('black')
        l.get_frame().set_alpha(None)

    # Formatting the limits, labels etc
    for a in [ax]:
        a.set_ylim([0, 9.5])
        a.set_xlim([0, 350])
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
        a.set_xticks([0, 100, 200, 300])
        a.set_xticklabels(('0', '100', '200', '300'))

plt.savefig('../../Figures/pdfs/elastic_grav_accs_codes.pdf')

plt.show()